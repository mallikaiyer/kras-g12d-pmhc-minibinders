"""Tests for the NVIDIA NIM fallback transport.

No network is touched: ``generate_one`` takes its opener as an argument
precisely so the retry, polling and error paths can be exercised offline. The
one piece of real logic here is ``relabel_chains``, and it is tested against a
synthetic single-chain array of exactly the shape RFdiffusion returns, because
every interface number computed downstream depends on it being right.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
from biotite.structure import AtomArray

from kras_pmhc import nim


class _Response:
    """Minimal stand-in for the object urlopen yields."""

    def __init__(self, status: int, body: dict, headers: dict | None = None):
        self.status = status
        self.headers = headers or {}
        self._body = json.dumps(body).encode()

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _opener(*responses):
    """An opener that returns the given responses in order."""
    queue = list(responses)
    calls = []

    def open_it(request):
        calls.append(request)
        return queue.pop(0)

    open_it.calls = calls
    return open_it


# --- key handling ---------------------------------------------------------


def test_a_missing_key_file_says_how_to_make_one(tmp_path):
    """With no file and no terminal, the error must name both routes."""
    with pytest.raises(nim.NimError, match="--set-key"):
        nim.read_api_key(tmp_path / "absent", allow_prompt=False)


def test_a_missing_key_is_prompted_for_when_there_is_a_terminal(tmp_path, monkeypatch):
    path = tmp_path / "key"
    monkeypatch.setattr(nim.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr(nim.getpass, "getpass", lambda _: "  nvapi-typed\n")

    assert nim.read_api_key(path) == "nvapi-typed"
    assert path.read_text().strip() == "nvapi-typed"


def test_a_prompted_key_is_not_world_readable(tmp_path, monkeypatch):
    """A key file must never be readable by anyone but its owner."""
    path = tmp_path / "key"
    monkeypatch.setattr(nim.getpass, "getpass", lambda _: "nvapi-typed")
    nim.prompt_for_key(path)
    assert path.stat().st_mode & 0o077 == 0


def test_a_mistyped_key_is_rejected_at_the_prompt(tmp_path, monkeypatch):
    path = tmp_path / "key"
    monkeypatch.setattr(nim.getpass, "getpass", lambda _: "sk-wrong-provider")
    with pytest.raises(nim.NimError, match="nvapi-"):
        nim.prompt_for_key(path)
    assert not path.exists(), "a rejected key must not be written"


def test_a_key_of_the_wrong_shape_is_rejected(tmp_path):
    path = tmp_path / "key"
    path.write_text("sk-not-an-nvidia-key")
    with pytest.raises(nim.NimError, match="nvapi-"):
        nim.read_api_key(path)


def test_a_good_key_is_stripped(tmp_path):
    path = tmp_path / "key"
    path.write_text("  nvapi-abc123\n")
    assert nim.read_api_key(path) == "nvapi-abc123"


# --- payload --------------------------------------------------------------


def test_only_coordinate_records_are_sent():
    pdb = "HEADER something\nATOM      1  N   ALA A   1\nEND\n"
    assert nim.atom_lines(pdb).splitlines() == ["ATOM      1  N   ALA A   1"]


def test_a_target_with_no_atoms_is_an_error():
    with pytest.raises(nim.NimError, match="no ATOM"):
        nim.atom_lines("HEADER only\n")


def test_the_payload_carries_the_agreed_parameters():
    payload = nim.build_payload(
        "ATOM      1  N   ALA A   1\n", "A1-180/0 C1-10/0 60-75",
        ["C5", "C6", "C7"], 50, 7,
    )
    assert payload["contigs"] == "A1-180/0 C1-10/0 60-75"
    assert payload["hotspot_res"] == ["C5", "C6", "C7"]
    assert payload["diffusion_steps"] == 50
    assert payload["random_seed"] == 7


def test_an_unseeded_payload_omits_the_field_entirely():
    payload = nim.build_payload("ATOM      1  N   ALA A   1\n", "A1-10/0 20", [], 50, None)
    assert "random_seed" not in payload


@pytest.mark.parametrize("steps", [0, 14, 51, 200])
def test_diffusion_steps_outside_the_documented_range_are_refused(steps):
    with pytest.raises(nim.NimError, match="documented"):
        nim.build_payload("ATOM      1  N   ALA A   1\n", "A1-10/0 20", [], steps, 0)


# --- HTTP paths -----------------------------------------------------------


def test_a_plain_success_returns_the_structure():
    opener = _opener(_Response(200, {"output_pdb": "ATOM...", "elapsed_ms": 1234}))
    design = nim.generate_one("nvapi-x", {"random_seed": 3}, opener)
    assert design.pdb_text == "ATOM..."
    assert design.seed == 3
    assert design.server_elapsed_ms == 1234


def test_a_202_is_polled_until_it_completes(monkeypatch):
    monkeypatch.setattr(nim, "POLL_INTERVAL_S", 0)
    opener = _opener(
        _Response(202, {}, {"nvcf-reqid": "abc"}),
        _Response(202, {}),
        _Response(200, {"output_pdb": "ATOM...", "elapsed_ms": 10}),
    )
    design = nim.generate_one("nvapi-x", {"random_seed": 0}, opener)
    assert design.pdb_text == "ATOM..."
    assert len(opener.calls) == 3


def test_a_202_without_a_request_id_cannot_be_polled():
    opener = _opener(_Response(202, {}))
    with pytest.raises(nim.NimError, match="nvcf-reqid"):
        nim.generate_one("nvapi-x", {}, opener)


def test_a_rate_limit_is_retried(monkeypatch):
    monkeypatch.setattr(nim, "BACKOFF_BASE_S", 0)
    opener = _opener(
        _Response(429, {"detail": "slow down"}),
        _Response(200, {"output_pdb": "ATOM...", "elapsed_ms": 5}),
    )
    design = nim.generate_one("nvapi-x", {}, opener)
    assert design.pdb_text == "ATOM..."


def test_a_bad_request_fails_immediately_rather_than_retrying():
    """A 400 means our contigs or target are wrong; retrying wastes the clock."""
    opener = _opener(_Response(400, {"detail": "bad contigs"}))
    with pytest.raises(nim.NimError, match="400"):
        nim.generate_one("nvapi-x", {}, opener)
    assert len(opener.calls) == 1


def test_a_success_without_a_structure_is_an_error():
    opener = _opener(_Response(200, {"elapsed_ms": 9}))
    with pytest.raises(nim.NimError, match="output_pdb"):
        nim.generate_one("nvapi-x", {}, opener)


# --- contigs and chain relabelling ---------------------------------------


def test_the_fixed_target_segments_are_read_from_the_contigs():
    assert nim.parse_contig_segments("A1-180/0 C1-10/0 60-75") == [
        ("A", 1, 180), ("C", 1, 10),
    ]


def test_a_contig_string_with_no_target_is_an_error():
    with pytest.raises(nim.NimError, match="no fixed target"):
        nim.parse_contig_segments("60-75")


def _single_chain(n_residues: int, atoms_per_residue: int = 4) -> AtomArray:
    """What RFdiffusion returns: one chain, renumbered from 1."""
    array = AtomArray(n_residues * atoms_per_residue)
    names = ["N", "CA", "C", "O"][:atoms_per_residue]
    for i in range(n_residues):
        for j, name in enumerate(names):
            k = i * atoms_per_residue + j
            array.chain_id[k] = "A"
            array.res_id[k] = i + 1
            array.res_name[k] = "GLY"
            array.atom_name[k] = name
            array.element[k] = name[0]
    array.coord = np.arange(len(array) * 3, dtype=float).reshape(-1, 3)
    return array


def test_the_target_keeps_its_chains_and_the_binder_gets_its_own():
    array = _single_chain(180 + 10 + 65)
    fixed = nim.relabel_chains(array, "A1-180/0 C1-10/0 60-75")

    hla = fixed[fixed.chain_id == "A"]
    peptide = fixed[fixed.chain_id == "C"]
    binder = fixed[fixed.chain_id == nim.BINDER_CHAIN_ID]

    assert len(np.unique(hla.res_id)) == 180
    assert len(np.unique(peptide.res_id)) == 10
    assert len(np.unique(binder.res_id)) == 65
    assert peptide.res_id.min() == 1 and peptide.res_id.max() == 10
    assert binder.res_id.min() == 1 and binder.res_id.max() == 65


def test_output_with_no_binder_is_rejected_rather_than_relabelled():
    """If the API read the contigs differently, we must not guess."""
    array = _single_chain(190)
    with pytest.raises(nim.NimError, match="no binder"):
        nim.relabel_chains(array, "A1-180/0 C1-10/0 60-75")


def test_relabelling_leaves_coordinates_untouched():
    array = _single_chain(180 + 10 + 60)
    fixed = nim.relabel_chains(array, "A1-180/0 C1-10/0 60-75")
    assert np.array_equal(fixed.coord, array.coord)
