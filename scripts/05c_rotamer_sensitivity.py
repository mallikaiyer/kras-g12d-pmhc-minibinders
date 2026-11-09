#!/usr/bin/env python
"""Stage 05 - how far the design gate depends on the constructed Asp6 rotamer.

7OW3 carries wild-type Gly at peptide position 6, so Asp6 does not exist in any
crystal structure of this complex at the conformer used here: stage 02 scanned
its rotamers, stage 03 built one, stage 04 minimised it. Every Asp6 number in
the stage-05 manifest is therefore measured against a side chain this project
constructed, and the gate's burial criterion partly measures that choice.

This script quantifies "partly". It re-scores every generated binder against an
independently built and relaxed alternative rotamer, changing nothing else, and
reports which gate verdicts move. It writes a table rather than altering any
decision: the project scores on the canonical rotamer, and this is the evidence
for the limitation that choice carries.

Read-only with respect to data/backbones and results/interface.

Usage:
    python scripts/05c_rotamer_sensitivity.py
    python scripts/05c_rotamer_sensitivity.py --alternative <other_crop.pdb>
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from biotite.structure.io.pdb import PDBFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import backbones, interface  # noqa: E402
from kras_pmhc.constants import MUTANT_POSITION  # noqa: E402

RELAXED = ROOT / "data" / "targets_relaxed"
MANIFEST = ROOT / "data" / "backbones" / "MANIFEST.tsv"
BACKBONE_DIR = ROOT / "data" / "backbones"
OUT = ROOT / "results" / "rotamer_sensitivity.tsv"

DEFAULT_CANONICAL = "7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb"
DEFAULT_ALTERNATIVE = "7ow3_L_G6D_chi1m180_chi2p0_relaxed_crop.pdb"

TARGET_CHAINS = ["A", "C"]
PEPTIDE_CHAIN = "C"
HLA_CHAIN = "A"

# The gate, repeated here deliberately rather than imported: this script must
# apply the same rule to the alternative rotamer that 05_ingest_backbones.py
# applied to the canonical one, and a change to either should be a visible
# change to both.
GATE_MIN_ASP6_SASA_LOSS = 10.0
GATE_MIN_PEPTIDE_SHARE = 0.10

FIELDS = (
    "backbone_id",
    "canonical_class", "canonical_sasa_loss", "canonical_share", "canonical_passes",
    "alternative_class", "alternative_sasa_loss", "alternative_share",
    "alternative_passes",
    "min_sasa_loss", "verdict_change",
)


def read(path: Path):
    return PDBFile.read(str(path)).get_structure(model=1)


def score_against(design, reference) -> tuple[str, float, float, bool]:
    """Apply the stage-05 gate to one design against one target rotamer.

    The design is superimposed by its *target* backbone, so swapping the
    reference does not depend on the two relaxed targets sharing a frame.
    """
    aligned, _ = backbones.superimpose_on_target(design, reference, TARGET_CHAINS)
    binder = backbones.binder_chain(aligned, TARGET_CHAINS)
    complexed = backbones.splice_binder(reference, aligned, binder)

    records = interface.classify(complexed, [binder], TARGET_CHAINS)
    areas = interface.interface_area_by_chain(records)
    peptide = areas.get(PEPTIDE_CHAIN, 0.0)
    total = peptide + areas.get(HLA_CHAIN, 0.0)
    share = peptide / total if total else 0.0

    asp6 = interface.find(records, PEPTIDE_CHAIN, MUTANT_POSITION)
    if asp6 is None:
        raise SystemExit(f"peptide residue {MUTANT_POSITION} absent from the complex")

    passes = (
        asp6.klass in interface.INTERFACE_CLASSES
        and asp6.sasa_loss >= GATE_MIN_ASP6_SASA_LOSS
        and share >= GATE_MIN_PEPTIDE_SHARE
    )
    return asp6.klass, asp6.sasa_loss, share, passes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--canonical", default=DEFAULT_CANONICAL)
    parser.add_argument("--alternative", default=DEFAULT_ALTERNATIVE)
    args = parser.parse_args()

    if not MANIFEST.is_file():
        raise SystemExit(f"{MANIFEST} not found - run scripts/05_ingest_backbones.py")

    alternative = RELAXED / args.alternative
    if not alternative.is_file():
        raise SystemExit(f"alternative target not found: {alternative}")

    rows = list(csv.DictReader(MANIFEST.open(), delimiter="\t"))
    reference = read(alternative)

    print(f"canonical  : {args.canonical}")
    print(f"alternative: {args.alternative}")
    print(f"re-scoring {len(rows)} design(s)\n")

    out_rows, flips = [], []
    for row in rows:
        design = read(BACKBONE_DIR / f"{row['backbone_id']}.pdb")
        klass, loss, share, passes = score_against(design, reference)

        was = row["passes_gate"] == "True"
        change = "" if was == passes else ("lost" if was else "gained")
        if change:
            flips.append((row["backbone_id"], change))

        out_rows.append({
            "backbone_id": row["backbone_id"],
            "canonical_class": row["asp6_class"],
            "canonical_sasa_loss": row["asp6_sasa_loss"],
            "canonical_share": row["peptide_share"],
            "canonical_passes": was,
            "alternative_class": klass,
            "alternative_sasa_loss": round(loss, 2),
            "alternative_share": round(share, 4),
            "alternative_passes": passes,
            "min_sasa_loss": round(min(float(row["asp6_sasa_loss"]), loss), 2),
            "verdict_change": change,
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    canonical_pass = {r["backbone_id"] for r in out_rows if r["canonical_passes"]}
    alternative_pass = {r["backbone_id"] for r in out_rows if r["alternative_passes"]}
    both = canonical_pass & alternative_pass

    for name, change in flips:
        print(f"  {name:<9} {change}")
    print(f"\npass set: {len(canonical_pass)} canonical, {len(alternative_pass)} "
          f"alternative, {len(both)} under both")
    print(f"verdicts moved: {len(flips)} of {len(rows)}")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    print("Selection uses the canonical rotamer; this table is the limitation's "
          "evidence, not a second gate.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
