"""NVIDIA NIM RFdiffusion client - the stage-05 fallback transport.

This module exists because Colab's free GPU quota ran out on 9 Nov 2026 with
the pilot still unrun. It is a *fallback*, and it is worse than the notebook in
one specific, recorded way: the hosted API exposes no noise-scale field, so
``denoiser.noise_scale_ca`` / ``noise_scale_frame`` cannot be set to the 0.5
that the RFdiffusion authors recommend for protein-protein interfaces. Designs
produced here are sampled at the default noise. That is a deviation, not a
preference; see docs/colab_runbook.md.

What the API does give us is scriptability: no browser session to babysit, no
GPU allocation to lose, and seeds we control. Everything downstream is
unchanged - the client writes the same directory layout the Colab notebook
wrote, so ``scripts/05_ingest_backbones.py`` reads it without modification and
the Levy (2010) interface classification is applied to these backbones exactly
as it would have been to Colab's.

Only the standard library is used for HTTP. One POST per design does not
justify a dependency, and keeping it stdlib means the script runs in any env
that can already import biotite.

Reference: https://docs.nvidia.com/nim/bionemo/rfdiffusion/latest/endpoints.html
"""

from __future__ import annotations

import getpass
import json
import os
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray

#: Hosted inference endpoint. The self-hosted NIM exposes the same path.
GENERATE_URL = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"
#: NVCF hands back a request id for long-running jobs; results are polled here.
STATUS_URL = "https://health.api.nvidia.com/v1/status/{req_id}"

#: Documented range for ``diffusion_steps``; 50 matches ``diffuser.T`` in the
#: notebook, so the two routes differ only in noise scale.
MIN_DIFFUSION_STEPS = 15
MAX_DIFFUSION_STEPS = 50

#: Retry policy. 429 and 5xx are transient on a free tier and worth waiting out;
#: 4xx other than 429 means the request itself is wrong and retrying is pointless.
RETRY_STATUSES = frozenset({429, 500, 502, 503, 504})
MAX_ATTEMPTS = 5
BACKOFF_BASE_S = 4.0

#: How long to wait between polls of a 202-accepted request, and how long to
#: keep polling before giving up on it.
POLL_INTERVAL_S = 5.0
POLL_TIMEOUT_S = 900.0

#: How closely the putative target CA trace must match the reference for a
#: chain ordering to be accepted. The correct ordering comes in well under
#: this (0.09 A observed); the wrong one is off by tens of angstroms.
LAYOUT_RMSD_TOLERANCE_A = 1.0

#: Chain id given to the generated binder. The target keeps A and C, matching
#: the staged crop, so the ingest script's chain assumptions still hold.
BINDER_CHAIN_ID = "B"


class NimError(RuntimeError):
    """A request failed in a way that is not worth retrying."""


@dataclass(frozen=True)
class Design:
    """One generated backbone and what produced it."""

    seed: int
    pdb_text: str
    elapsed_s: float
    server_elapsed_ms: int | None


#: Every key issued by build.nvidia.com carries this prefix.
KEY_PREFIX = "nvapi-"


def _validate_key(key: str, source: str) -> str:
    if not key.startswith(KEY_PREFIX):
        raise NimError(f"{source} does not look like an NVIDIA key (expected nvapi-...)")
    return key


def prompt_for_key(path: Path) -> str:
    """Ask for the key at the terminal and save it, readable only by its owner.

    ``getpass`` is used rather than ``input`` so the key is never echoed, and
    the file is created with mode 0600 from the outset rather than chmodded
    afterwards - a key must not be briefly world-readable.
    """
    key = _validate_key(getpass.getpass("NVIDIA API key (nvapi-...): ").strip(), "input")
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with open(handle, "w") as out:
        out.write(key + "\n")
    print(f"saved to {path} (mode 600, gitignored)")
    return key


def read_api_key(path: Path, allow_prompt: bool = True) -> str:
    """Read the key from a file, or ask for it if there is a terminal to ask at.

    A key passed as a command-line argument lands in the shell history and in
    the process table; a key in a gitignored file or typed at a hidden prompt
    does neither.
    """
    if path.is_file():
        return _validate_key(path.read_text().strip(), str(path))
    if allow_prompt and sys.stdin.isatty():
        return prompt_for_key(path)
    raise NimError(
        f"no API key at {path}, and no terminal to ask at. Either run\n"
        f"  python scripts/05b_nim_generate.py --set-key\n"
        f"from a terminal, or copy the key from build.nvidia.com and run\n"
        f"  pbpaste > {path} && chmod 600 {path}"
    )


def atom_lines(pdb_text: str) -> str:
    """Keep only coordinate records.

    NVIDIA's own examples post ATOM lines alone. Headers and CRYST1 records are
    not needed to define the target and only inflate a request that has a size
    limit we do not know.
    """
    kept = [line for line in pdb_text.splitlines() if line.startswith(("ATOM", "TER"))]
    if not kept:
        raise NimError("target PDB contains no ATOM records")
    return "\n".join(kept) + "\n"


def build_payload(
    pdb_text: str,
    contigs: str,
    hotspots: list[str],
    diffusion_steps: int,
    seed: int | None,
) -> dict:
    """Assemble one request.

    ``hotspot_res`` is sent as a JSON array. The endpoint reference shows a
    comma-joined string in one place and an array in another; the quickstart
    uses an array, and a list is what the field's declared type says.
    """
    if not MIN_DIFFUSION_STEPS <= diffusion_steps <= MAX_DIFFUSION_STEPS:
        raise NimError(
            f"diffusion_steps={diffusion_steps} outside the documented "
            f"{MIN_DIFFUSION_STEPS}-{MAX_DIFFUSION_STEPS}"
        )
    payload = {
        "input_pdb": atom_lines(pdb_text),
        "contigs": contigs,
        "hotspot_res": hotspots,
        "diffusion_steps": diffusion_steps,
    }
    if seed is not None:
        payload["random_seed"] = seed
    return payload


def _request(url: str, key: str, payload: dict | None, opener) -> tuple[int, dict, dict]:
    """One HTTP call. Returns (status, headers, parsed body).

    ``opener`` is injected so tests exercise this module without a network.
    """
    headers = {
        "Authorization": f"Bearer {key}",
        "Accept": "application/json",
    }
    data = None
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"

    request = urllib.request.Request(url, data=data, headers=headers)
    try:
        with opener(request) as response:
            body = response.read().decode()
            parsed = json.loads(body) if body.strip() else {}
            return response.status, dict(response.headers), parsed
    except urllib.error.HTTPError as error:
        body = error.read().decode(errors="replace")
        try:
            parsed = json.loads(body) if body.strip() else {}
        except json.JSONDecodeError:
            parsed = {"raw": body[:2000]}
        return error.code, dict(error.headers or {}), parsed


def _poll(req_id: str, key: str, opener) -> dict:
    """Wait out a 202-accepted request.

    NVCF accepts long jobs and hands back a request id. Polling is capped so a
    job lost on the server side cannot hang the run forever.
    """
    deadline = time.monotonic() + POLL_TIMEOUT_S
    while time.monotonic() < deadline:
        time.sleep(POLL_INTERVAL_S)
        status, _, body = _request(STATUS_URL.format(req_id=req_id), key, None, opener)
        if status == 200:
            return body
        if status != 202:
            raise NimError(f"polling {req_id} returned {status}: {body}")
    raise NimError(f"request {req_id} did not finish within {POLL_TIMEOUT_S:.0f}s")


def generate_one(
    key: str,
    payload: dict,
    opener=urllib.request.urlopen,
) -> Design:
    """Post one design, retrying transient failures.

    The seed is read back out of the payload rather than passed separately, so
    the logged seed is necessarily the one that was sent.
    """
    seed = payload.get("random_seed")
    started = time.monotonic()

    for attempt in range(1, MAX_ATTEMPTS + 1):
        status, headers, body = _request(GENERATE_URL, key, payload, opener)

        if status == 202:
            req_id = headers.get("nvcf-reqid") or headers.get("NVCF-REQID")
            if not req_id:
                raise NimError("202 accepted but no nvcf-reqid header to poll")
            body = _poll(req_id, key, opener)
            status = 200

        if status == 200:
            pdb_text = body.get("output_pdb")
            if not pdb_text:
                raise NimError(f"200 OK but no output_pdb in response: {list(body)}")
            return Design(
                seed=seed,
                pdb_text=pdb_text,
                elapsed_s=time.monotonic() - started,
                server_elapsed_ms=body.get("elapsed_ms"),
            )

        if status in RETRY_STATUSES and attempt < MAX_ATTEMPTS:
            wait = BACKOFF_BASE_S * (2 ** (attempt - 1))
            print(f"    HTTP {status}; retrying in {wait:.0f}s (attempt {attempt})")
            time.sleep(wait)
            continue

        raise NimError(f"HTTP {status}: {json.dumps(body)[:1000]}")

    raise NimError(f"gave up after {MAX_ATTEMPTS} attempts")


def parse_contig_segments(contigs: str) -> list[tuple[str, int, int]]:
    """Split a contig string into its fixed target segments.

    ``"A1-180/0 C1-10/0 60-75"`` gives ``[("A", 1, 180), ("C", 1, 10)]``. The
    trailing binder range is deliberately not returned: its length is sampled
    per design and is only known from the output.
    """
    segments = []
    for token in contigs.split():
        token = token.removesuffix("/0")
        if not token or token[0].isdigit():
            continue  # the binder length range
        chain, span = token[0], token[1:]
        low, _, high = span.partition("-")
        segments.append((chain, int(low), int(high)))
    if not segments:
        raise NimError(f"no fixed target segments in contigs {contigs!r}")
    return segments


def _residue_starts(array: AtomArray) -> np.ndarray:
    """Index of the first atom of each residue, plus a terminating index."""
    boundaries = (
        np.where(
            (array.res_id[1:] != array.res_id[:-1])
            | (array.chain_id[1:] != array.chain_id[:-1])
        )[0]
        + 1
    )
    return np.concatenate(([0], boundaries, [len(array)]))


def resolve_layout(array: AtomArray, reference: AtomArray, contigs: str) -> bool:
    """Work out whether the generated chain precedes the target in the output.

    The contig string puts the binder last, but this endpoint's RFdiffusion
    writes the diffused chain **first**, and nothing in the response says so.
    Guessing is not acceptable here: labelling the binder's residues as HLA
    produced a plausible-looking manifest whose interface numbers were entirely
    fictitious, which is how this function came to exist.

    So the layout is measured. Each hypothesis is tested by superimposing the
    putative target CA trace onto the reference: the correct one reproduces the
    target's own geometry, because the target is held rigid during denoising.
    Observed on 9 Nov 2026: 0.09 A for binder-first against 21.70 A for
    target-first.

    Returns:
        True if the generated chain comes first.

    Raises:
        NimError: if neither hypothesis recovers the target geometry, which
            means the output does not contain the target we conditioned on.
    """
    segments = parse_contig_segments(contigs)
    chains = [chain for chain, _, _ in segments]
    n_target = sum(high - low + 1 for _, low, high in segments)

    generated = array[array.atom_name == "CA"]
    target = reference[
        (reference.atom_name == "CA") & np.isin(reference.chain_id, chains)
    ]
    if len(target) != n_target:
        raise NimError(
            f"reference has {len(target)} CA atoms in chains {chains} but the "
            f"contigs describe {n_target}"
        )
    n_binder = len(generated) - n_target
    if n_binder <= 0:
        raise NimError(
            f"output has {len(generated)} residues but the target alone needs "
            f"{n_target}; no binder was generated"
        )

    scores = {}
    for binder_first, window in (
        (False, slice(0, n_target)),
        (True, slice(n_binder, n_binder + n_target)),
    ):
        fitted, _ = struc.superimpose(target, generated[window])
        scores[binder_first] = float(struc.rmsd(target, fitted))

    binder_first = min(scores, key=scores.get)
    if scores[binder_first] > LAYOUT_RMSD_TOLERANCE_A:
        raise NimError(
            "neither chain ordering recovers the target geometry "
            f"(binder-first {scores[True]:.2f} A, target-first "
            f"{scores[False]:.2f} A); the output does not contain the target "
            "that was sent"
        )
    return binder_first


def relabel_chains(
    array: AtomArray, contigs: str, binder_first: bool = False
) -> AtomArray:
    """Restore the target's chain ids and numbering in RFdiffusion output.

    RFdiffusion renumbers its output as one continuous chain, which is why the
    Colab route called ColabDesign's ``fix_pdb``. This is the same operation
    done locally.

    ``binder_first`` says which end of the output the generated chain sits at.
    Do not set it by assumption - :func:`resolve_layout` measures it against
    the reference structure.

    Fails rather than guesses if the residue count does not match the contigs.
    A mismatch means the API interpreted the contig string differently than we
    did, and silently mislabelling chains would corrupt every interface number
    computed afterwards.
    """
    segments = parse_contig_segments(contigs)
    expected_target = sum(high - low + 1 for _, low, high in segments)

    starts = _residue_starts(array)
    n_residues = len(starts) - 1

    if n_residues <= expected_target:
        raise NimError(
            f"output has {n_residues} residues but the target alone needs "
            f"{expected_target}; no binder was generated"
        )

    target_labels = [
        (chain, res_id)
        for chain, low, high in segments
        for res_id in range(low, high + 1)
    ]
    binder_labels = [
        (BINDER_CHAIN_ID, i + 1) for i in range(n_residues - expected_target)
    ]
    labels = (
        binder_labels + target_labels if binder_first else target_labels + binder_labels
    )

    relabelled = array.copy()
    for index, (chain, res_id) in enumerate(labels):
        lo, hi = starts[index], starts[index + 1]
        relabelled.chain_id[lo:hi] = chain
        relabelled.res_id[lo:hi] = res_id
    return relabelled
