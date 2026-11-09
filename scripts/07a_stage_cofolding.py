#!/usr/bin/env python
"""Stage 07a - build Boltz-2 cofolding inputs for the selected sequences.

Three predictions per sequence, which are the three questions the project's
ranking asks, in order:

    g12d   binder + HLA alpha1/alpha2 + VVVGADGVGK   does it bind the target?
    wt     binder + HLA alpha1/alpha2 + VVVGAGGVGK   does it bind without Asp6?
    apo    binder + HLA alpha1/alpha2                does it bind the empty groove?

Mutation specificity is g12d against wt - the peptides differ at position 6 only.
Peptide-dependent binding is g12d against apo. Affinity is last, and comes from
the g12d prediction alone. A binder that scores equally on all three is an HLA
binder, not a neoepitope binder, however tightly it binds.

Cofolding is what makes this measurable at all: PRODIGY and every other contact
based method needs full-atom side chains on both partners, and stage 05 produces
backbones only. Predicting the complex from sequence supplies them.

The HLA chain's MSA is computed once in the notebook and referenced by every
input, rather than queried per prediction: the sequence is identical across all
of them, and the public MSA server should not be asked the same question 40
times. Each YAML carries the token __HLA_MSA__ where that path goes, and the
notebook substitutes it.

Usage:
    python scripts/07a_stage_cofolding.py --drive-dir "<...>/kras_pmhc_stage07"
    python scripts/07a_stage_cofolding.py --drive-dir "<...>" --per-design 2
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc.constants import (  # noqa: E402
    HLA_A1101_HEAVY_CHAIN,
    HLA_DESIGN_TARGET_RESIDUES,
)

SELECTED = ROOT / "data" / "sequences" / "SELECTED.tsv"

#: The design target's groove: heavy-chain residues 1-180, exactly the chain A
#: that every stage from 03 onwards has used. Verified against the structure by
#: tests/test_constants.py, so this slice is not a guess.
HLA_GROOVE = HLA_A1101_HEAVY_CHAIN[:HLA_DESIGN_TARGET_RESIDUES]

#: G12D decamer (KRAS 7-16) and its wild type. They differ at position 6 alone,
#: which is what makes the pair a specificity test rather than two experiments.
PEPTIDE_G12D = "VVVGADGVGK"
PEPTIDE_WT = "VVVGAGGVGK"

#: Chain ids match the design target: A groove, C peptide, B binder.
CHAIN_HLA, CHAIN_PEPTIDE, CHAIN_BINDER = "A", "C", "B"

#: Substituted by the notebook once the HLA MSA has been computed.
MSA_TOKEN = "__HLA_MSA__"

#: De novo binders have no homologues and a 10-mer has no useful alignment, so
#: both run single-sequence. Only the HLA gets an MSA.
MSA_EMPTY = "empty"

CONDITIONS = ("g12d", "wt", "apo")

MANIFEST_FIELDS = (
    "job_id", "backbone_id", "sample", "condition", "binder_sequence",
    "peptide_sequence", "mpnn_score", "yaml_file", "sha256", "staged_utc",
)


def load_selected(per_design: int) -> list[dict]:
    """The best sequences per design, by ProteinMPNN score (lower is better)."""
    if not SELECTED.is_file():
        raise SystemExit(
            f"{SELECTED} not found - run scripts/06b_select_sequences.py"
        )
    with SELECTED.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    by_design = collections.defaultdict(list)
    for row in rows:
        by_design[row["backbone_id"]].append(row)

    chosen = []
    for design in sorted(by_design):
        ranked = sorted(by_design[design], key=lambda r: float(r["score"]))
        chosen.extend(ranked[:per_design])
    return chosen


def yaml_for(binder: str, condition: str) -> str:
    """One Boltz-2 input.

    Written as text rather than through a YAML library deliberately: the file is
    small, its shape is fixed, and a reader auditing what was predicted should
    see exactly the bytes that were sent.
    """
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        f"      id: {CHAIN_HLA}",
        f"      sequence: {HLA_GROOVE}",
        f"      msa: {MSA_TOKEN}",
    ]
    if condition != "apo":
        peptide = PEPTIDE_G12D if condition == "g12d" else PEPTIDE_WT
        lines += [
            "  - protein:",
            f"      id: {CHAIN_PEPTIDE}",
            f"      sequence: {peptide}",
            f"      msa: {MSA_EMPTY}",
        ]
    lines += [
        "  - protein:",
        f"      id: {CHAIN_BINDER}",
        f"      sequence: {binder}",
        f"      msa: {MSA_EMPTY}",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--drive-dir", required=True, type=Path)
    parser.add_argument(
        "--per-design",
        type=int,
        default=1,
        help="best N sequences per design, by MPNN score (default 1). Each one "
             "costs three predictions.",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    chosen = load_selected(args.per_design)
    drive = args.drive_dir.expanduser()
    inputs = drive / "inputs"
    if inputs.exists() and not args.force:
        existing = list(inputs.glob("*.yaml"))
        if existing:
            raise SystemExit(
                f"{inputs} already holds {len(existing)} input(s); pass --force "
                "to replace them"
            )
    if args.force and inputs.exists():
        shutil.rmtree(inputs)
    inputs.mkdir(parents=True, exist_ok=True)

    rows = []
    for row in chosen:
        for condition in CONDITIONS:
            job_id = f"{row['backbone_id']}_s{row['sample']}_{condition}"
            text = yaml_for(row["sequence"], condition)
            path = inputs / f"{job_id}.yaml"
            path.write_text(text)
            rows.append({
                "job_id": job_id,
                "backbone_id": row["backbone_id"],
                "sample": row["sample"],
                "condition": condition,
                "binder_sequence": row["sequence"],
                "peptide_sequence": (
                    "" if condition == "apo"
                    else PEPTIDE_G12D if condition == "g12d" else PEPTIDE_WT
                ),
                "mpnn_score": row["score"],
                "yaml_file": path.name,
                "sha256": hashlib.sha256(text.encode()).hexdigest(),
                "staged_utc": datetime.now(UTC).isoformat(timespec="seconds"),
            })

    manifest = drive / "JOBS.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Also keep a copy in the repository: the Drive folder is scratch space, and
    # which sequences were predicted is part of the record.
    local = ROOT / "data" / "cofolding"
    local.mkdir(parents=True, exist_ok=True)
    shutil.copy2(manifest, local / "JOBS.tsv")

    designs = {r["backbone_id"] for r in rows}
    print(f"{len(chosen)} sequence(s) from {len(designs)} design(s)")
    print(f"{len(rows)} predictions: {len(chosen)} x {len(CONDITIONS)} conditions")
    print(f"HLA groove: {len(HLA_GROOVE)} residues")
    print(f"\nwrote {inputs}")
    print(f"wrote {manifest}")
    print(f"wrote {(local / 'JOBS.tsv').relative_to(ROOT)}")
    print("\nLet Drive finish syncing, then run notebooks/07_boltz2_cofolding.ipynb")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
