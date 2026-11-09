"""Tests for ProteinMPNN sequence design.

Offline: ``design`` takes its opener as an argument so the HTTP path is
exercised without a network. The multi-FASTA fixture below is copied from a
real response, which matters - the leading ``>input`` record and the exact
header format are the two things a hand-written fixture would get wrong.
"""

from __future__ import annotations

import json

import pytest

from kras_pmhc import mpnn

#: Verbatim from a response for nim3_6 on 9 Nov 2026, truncated to two designs.
# ruff: noqa: E501 - the header is real output; wrapping it would make the
# fixture unfaithful, and its exact format is what these tests check.
REAL_RESPONSE = """>input, score=1.7525, global_score=1.5131, fixed_chains=['A', 'C'], designed_chains=['B'], model_name=v_48_002, git_hash=unknown, seed=674
GGGGGGGGGG
>T=0.1, sample=1, score=1.1496, global_score=1.3254, seq_recovery=0.0143
SLEEEEEARE
>T=0.1, sample=2, score=1.1970, global_score=1.3362, seq_recovery=0.0143
SLLAALQAQK
"""

PDB = "ATOM      1  N   GLY B   1      0.000   0.000   0.000\n"


class _Response:
    def __init__(self, status, body):
        self.status = status
        self.headers = {}
        self._body = json.dumps(body).encode()

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _opener(body, status=200):
    def open_it(request):
        return _Response(status, body)

    return open_it


# --- parsing --------------------------------------------------------------


def test_the_input_record_is_not_counted_as_a_design():
    """The first record is the backbone's own sequence - all glycine here.

    Counting it would put a poly-glycine 'design' into the submission.
    """
    designs = mpnn.parse_mfasta(REAL_RESPONSE)
    assert len(designs) == 2
    assert all("G" * 10 != d.sequence for d in designs)


def test_the_header_fields_are_read_back():
    first = mpnn.parse_mfasta(REAL_RESPONSE)[0]
    assert first.sample == 1
    assert first.temperature == 0.1
    assert first.score == pytest.approx(1.1496)
    assert first.global_score == pytest.approx(1.3254)
    assert first.seq_recovery == pytest.approx(0.0143)
    assert first.sequence == "SLEEEEEARE"
    assert first.length == 10


def test_a_wrong_length_is_an_error_not_a_warning():
    """A length mismatch means the wrong chain was designed.

    Accepting it would attach a sequence to the wrong residues, and every
    downstream claim would be about a protein that was never designed.
    """
    with pytest.raises(mpnn.MpnnError, match="wrong chain"):
        mpnn.parse_mfasta(REAL_RESPONSE, expected_length=70)


def test_the_expected_length_passes_when_it_matches():
    designs = mpnn.parse_mfasta(REAL_RESPONSE, expected_length=10)
    assert len(designs) == 2


def test_a_response_with_only_the_input_record_is_an_error():
    only_input = "\n".join(REAL_RESPONSE.splitlines()[:2])
    with pytest.raises(mpnn.MpnnError, match="no designed sequences"):
        mpnn.parse_mfasta(only_input)


def test_a_header_with_no_sequence_is_an_error():
    truncated = REAL_RESPONSE + ">T=0.1, sample=3, score=1.0\n"
    with pytest.raises(mpnn.MpnnError, match="no sequence"):
        mpnn.parse_mfasta(truncated)


# --- payload --------------------------------------------------------------


def test_the_payload_designs_only_the_named_chain():
    payload = mpnn.build_payload(PDB, ["B"], num_sequences=10)
    assert payload["input_pdb_chains"] == ["B"]
    assert payload["num_seq_per_target"] == 10
    assert payload["sampling_temp"] == [0.1]


def test_the_soluble_model_is_off():
    """Measured, not assumed: it gave 71-80% hydrophobic sequences here."""
    assert mpnn.build_payload(PDB, ["B"])["use_soluble_model"] is False


def test_designing_no_chain_is_an_error():
    with pytest.raises(mpnn.MpnnError, match="no chain"):
        mpnn.build_payload(PDB, [])


@pytest.mark.parametrize("n", [0, -1, 101])
def test_sequence_counts_outside_the_documented_range_are_refused(n):
    with pytest.raises(mpnn.MpnnError, match="1-100"):
        mpnn.build_payload(PDB, ["B"], num_sequences=n)


@pytest.mark.parametrize("temperature", [0.0, -0.1, 1.5])
def test_temperatures_outside_the_documented_range_are_refused(temperature):
    with pytest.raises(mpnn.MpnnError, match=r"\(0, 1\]"):
        mpnn.build_payload(PDB, ["B"], temperature=temperature)


def test_a_seed_is_sent_when_given_even_though_it_is_ignored():
    """Recorded as a claim about the request, not about reproducibility."""
    assert mpnn.build_payload(PDB, ["B"], seed=7)["random_seed"] == 7
    assert "random_seed" not in mpnn.build_payload(PDB, ["B"])


# --- the call -------------------------------------------------------------


def test_design_returns_the_raw_response_alongside_the_parse():
    """The caller must be able to keep the raw text: the seed is ignored, so
    the response cannot be regenerated."""
    designs, raw = mpnn.design(
        "nvapi-x", PDB, ["B"], expected_length=10, opener=_opener({"mfasta": REAL_RESPONSE})
    )
    assert len(designs) == 2
    assert raw == REAL_RESPONSE


def test_a_response_without_sequences_is_an_error():
    with pytest.raises(mpnn.MpnnError, match="no mfasta"):
        mpnn.design("nvapi-x", PDB, ["B"], opener=_opener({"scores": [1.0]}))
