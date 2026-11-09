#!/usr/bin/env python
"""Stage 07b - cofold the selected designs via the hosted Boltz-2 endpoint.

A backup for notebooks/07_boltz2_cofolding.ipynb, run in parallel with it rather
than instead of it. The two differ in one scientifically relevant way:

* The Colab run supplies an MSA for the HLA groove.
* This one runs **single-sequence** for every chain, because using an alignment
  here would mean calling the MSA-search endpoint as well.

Single-sequence prediction is weaker for the HLA, which is a well-aligned
protein with thousands of homologues. It is adequate for a backup and makes the
two runs independent: agreement between them is then evidence, not shared
method. Any number taken from this run must say which run it came from.

Scope, set by the user: the best sequence per design by ProteinMPNN score, each
folded against the G12D and wild-type peptides. No apo condition.

    g12d   HLA alpha1/alpha2 + VVVGADGVGK + binder
    wt     HLA alpha1/alpha2 + VVVGAGGVGK + binder

The pair differs at peptide position 6 alone, so the comparison isolates the
mutation - which is the project's first ranking criterion.

Writes after every prediction, so an interruption costs one job:

    data/cofolding_nim/<job_id>/model.cif
    data/cofolding_nim/<job_id>/confidence.json
    data/cofolding_nim/RESULTS.tsv

Usage:
    python scripts/07b_nim_cofolding.py
    python scripts/07b_nim_cofolding.py --per-design 2
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
import time
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import nim  # noqa: E402
from kras_pmhc.constants import (  # noqa: E402
    HLA_A1101_HEAVY_CHAIN,
    HLA_DESIGN_TARGET_RESIDUES,
)

SELECTED = ROOT / "data" / "sequences" / "SELECTED.tsv"
OUT_DIR = ROOT / "data" / "cofolding_nim"
RESULTS = OUT_DIR / "RESULTS.tsv"
KEY_FILE = ROOT / ".nim_api_key"

HLA_GROOVE = HLA_A1101_HEAVY_CHAIN[:HLA_DESIGN_TARGET_RESIDUES]
PEPTIDE = {"g12d": "VVVGADGVGK", "wt": "VVVGAGGVGK"}
CONDITIONS = ("g12d", "wt")

#: Chain ids match the design target throughout: A groove, C peptide, B binder.
CHAINS = {"hla": "A", "peptide": "C", "binder": "B"}

#: Endpoint defaults, stated explicitly so the printed record is complete.
RECYCLING_STEPS = 3
SAMPLING_STEPS = 50
DIFFUSION_SAMPLES = 1
STEP_SCALE = 1.638

FIELDS = (
    "job_id", "backbone_id", "sample", "condition", "binder_sequence",
    "peptide_sequence", "mpnn_score", "iptm", "ptm", "protein_iptm",
    "complex_plddt", "complex_iplddt", "confidence_score",
    "pair_chains_iptm", "elapsed_s", "predicted_utc",
)


def best_per_design(per_design: int) -> list[dict]:
    if not SELECTED.is_file():
        raise SystemExit(f"{SELECTED} not found - run scripts/06b_select_sequences.py")
    with SELECTED.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    grouped = collections.defaultdict(list)
    for row in rows:
        grouped[row["backbone_id"]].append(row)

    chosen = []
    for design in sorted(grouped):
        chosen.extend(
            sorted(grouped[design], key=lambda r: float(r["score"]))[:per_design]
        )
    return chosen


def payload_for(binder: str, condition: str) -> dict:
    """One request. No msa field anywhere, so every chain is single-sequence."""
    polymers = [
        {
            "id": CHAINS["hla"],
            "molecule_type": "protein",
            "sequence": HLA_GROOVE,
        },
        {
            "id": CHAINS["peptide"],
            "molecule_type": "protein",
            "sequence": PEPTIDE[condition],
        },
        {
            "id": CHAINS["binder"],
            "molecule_type": "protein",
            "sequence": binder,
        },
    ]
    return {
        "polymers": polymers,
        "recycling_steps": RECYCLING_STEPS,
        "sampling_steps": SAMPLING_STEPS,
        "diffusion_samples": DIFFUSION_SAMPLES,
        "step_scale": STEP_SCALE,
        "without_potentials": False,
        "output_format": "mmcif",
    }


def first(values, default=""):
    """The top-ranked diffusion sample's score, where the API returns a list."""
    if isinstance(values, list) and values:
        return values[0]
    return values if values is not None else default


def append(row: dict) -> None:
    exists = RESULTS.is_file()
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    with RESULTS.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def already_done() -> set[str]:
    if not RESULTS.is_file():
        return set()
    with RESULTS.open() as handle:
        return {row["job_id"] for row in csv.DictReader(handle, delimiter="\t")}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--per-design", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    key = nim.read_api_key(KEY_FILE, allow_prompt=False)
    chosen = best_per_design(args.per_design)
    done = set() if args.force else already_done()

    jobs = [
        (row, condition)
        for row in chosen
        for condition in CONDITIONS
        if f"{row['backbone_id']}_s{row['sample']}_{condition}" not in done
    ]

    print(f"{len(chosen)} sequence(s) x {len(CONDITIONS)} conditions "
          f"= {len(chosen) * len(CONDITIONS)} predictions")
    print("single-sequence mode (no MSA) - this is the backup run")
    print(f"{len(done)} already done, {len(jobs)} to run\n")

    failures = []
    for index, (row, condition) in enumerate(jobs, start=1):
        job_id = f"{row['backbone_id']}_s{row['sample']}_{condition}"
        print(f"[{index}/{len(jobs)}] {job_id}", flush=True)

        started = time.monotonic()
        try:
            body = nim.post(nim.BOLTZ2_URL, key, payload_for(row["sequence"], condition))
        except nim.NimError as error:
            print(f"  FAILED: {error}", flush=True)
            failures.append(job_id)
            continue
        elapsed = time.monotonic() - started

        structures = body.get("structures") or []
        if not structures:
            print("  FAILED: no structure in response", flush=True)
            failures.append(job_id)
            continue

        out = OUT_DIR / job_id
        out.mkdir(parents=True, exist_ok=True)
        (out / "model.cif").write_text(structures[0]["structure"])
        # The whole response minus the structures: every score the service
        # returned, kept verbatim rather than only the fields tabulated below.
        (out / "confidence.json").write_text(
            json.dumps({k: v for k, v in body.items() if k != "structures"}, indent=1)
        )

        append({
            "job_id": job_id,
            "backbone_id": row["backbone_id"],
            "sample": row["sample"],
            "condition": condition,
            "binder_sequence": row["sequence"],
            "peptide_sequence": PEPTIDE[condition],
            "mpnn_score": row["score"],
            "iptm": first(body.get("iptm_scores")),
            "ptm": first(body.get("ptm_scores")),
            "protein_iptm": first(body.get("protein_iptm_scores")),
            "complex_plddt": first(body.get("complex_plddt_scores")),
            "complex_iplddt": first(body.get("complex_iplddt_scores")),
            "confidence_score": first(body.get("confidence_scores")),
            "pair_chains_iptm": json.dumps(first(body.get("pair_chains_iptm_scores"), {})),
            "elapsed_s": round(elapsed, 1),
            "predicted_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        })
        print(f"  {elapsed:.0f}s  iptm {first(body.get('iptm_scores'))}  "
              f"plddt {first(body.get('complex_plddt_scores'))}", flush=True)

    print(f"\n{len(already_done())}/{len(chosen) * len(CONDITIONS)} complete")
    if failures:
        print(f"failed: {failures}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
