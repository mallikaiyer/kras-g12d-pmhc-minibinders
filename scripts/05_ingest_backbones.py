#!/usr/bin/env python
"""Stage 05 - ingest RFdiffusion backbones from Drive and apply the pilot gate.

Reads the backbones generated in Colab straight from the synced Google Drive
folder, copies them into the repository, and scores each one. The question this
stage answers is not "is this a binder" - nothing here can know that - but:

  * did RFdiffusion leave the target coordinates alone,
  * is the generated chain the length we asked for, and
  * does the binder actually engage Asp6, or only the alpha1/alpha2 helix tops?

That last question is the whole point of the hotspots, and it is the predicted
failure mode: the helix tops are a far larger hydrophobic surface than a
ten-residue peptide lying in a groove.

Interfaces are identified by the Levy (2010) rSASA classification in
``kras_pmhc.interface``, not by a distance cutoff. The closest-approach distance
to the carboxylate is reported for description only and decides nothing.

Usage:
    python scripts/05_ingest_backbones.py --drive-dir "<...>/My Drive/kras_pmhc_stage05"
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import backbones, interface, nim, structures  # noqa: E402
from kras_pmhc.constants import ASP_CARBOXYLATE_ATOMS, MUTANT_POSITION  # noqa: E402
from kras_pmhc.geometry import min_distance  # noqa: E402

RELAXED = ROOT / "data" / "targets_relaxed"
OUT_DIR = ROOT / "data" / "backbones"
INTERFACE_DIR = ROOT / "results" / "interface"

#: Chains of the staged target. Chain A is HLA alpha1/alpha2, C the peptide.
TARGET_CHAINS = ["A", "C"]
PEPTIDE_CHAIN = "C"
HLA_CHAIN = "A"

# --- pilot gate -----------------------------------------------------------
# Both conditions must hold for a backbone to count as engaging the mutation.
#: Asp6 must be an interface residue - support, rim or core.
GATE_ASP6_CLASSES = interface.INTERFACE_CLASSES
#: and the peptide must carry at least this share of target-side interface area.
GATE_MIN_PEPTIDE_SHARE = 0.20
#: Proceed to production only if at least this many of the pilot's 8 pass.
GATE_MIN_PASSING = 4

# --- hard validation ------------------------------------------------------
#: RFdiffusion holds the target fixed, so this is a parsing check, not physics.
MAX_TARGET_RMSD = 0.1
MIN_BINDER_LENGTH = 60
MAX_BINDER_LENGTH = 75
#: Below this the generated chain overlaps the target rather than touching it.
#: Measured against the *whole* target, not the carboxylate: a binder can sit a
#: comfortable distance from Asp6 while interpenetrating the HLA elsewhere, and
#: that failure is invisible to a carboxylate-only check.
MIN_HEAVY_ATOM_APPROACH = 2.0

MANIFEST_FIELDS = (
    "backbone_id", "batch", "seed", "binder_chain", "binder_length",
    "target_rmsd", "asp6_class", "asp6_delta_rsasa", "asp6_sasa_loss",
    "peptide_interface_area", "hla_interface_area", "peptide_share",
    "n_core_binder", "n_support_binder", "n_rim_binder",
    "n_interface_binder", "n_interface_target", "buried_interface_area",
    "helix", "strand", "coil",
    "min_dist_to_target", "min_dist_to_carboxylate", "passes_gate",
)


def read_pdb(path: Path):
    return PDBFile.read(str(path)).get_structure(model=1)


def load_reference(target_name: str):
    """The staged target, as the generated structures should still contain it."""
    path = RELAXED / target_name
    if not path.is_file():
        raise SystemExit(f"reference target not found: {path}")
    return read_pdb(path)


def read_run_log(drive_dir: Path) -> dict[str, dict]:
    path = drive_dir / "run_log.tsv"
    if not path.is_file():
        raise SystemExit(
            f"{path} not found - has the notebook run, and has Drive finished "
            "syncing back?"
        )
    with path.open() as handle:
        return {row["filename"]: row for row in csv.DictReader(handle, delimiter="\t")}


def score(array, reference, backbone_id: str, log_row: dict) -> tuple[dict, list]:
    """Metrics for one generated backbone, plus its per-residue Levy table."""
    binder = backbones.binder_chain(array, TARGET_CHAINS)

    records = interface.classify(array, [binder], TARGET_CHAINS)
    areas = interface.interface_area_by_chain(records)
    peptide_area = areas.get(PEPTIDE_CHAIN, 0.0)
    hla_area = areas.get(HLA_CHAIN, 0.0)
    target_area = peptide_area + hla_area

    asp6 = interface.find(records, PEPTIDE_CHAIN, MUTANT_POSITION)
    if asp6 is None:
        raise ValueError(
            f"{backbone_id}: peptide residue {MUTANT_POSITION} absent; the "
            "target was not carried through generation intact"
        )

    # Share of the *target* side only: the binder's own buried area says
    # nothing about which part of the target it chose.
    peptide_share = peptide_area / target_area if target_area > 0 else 0.0

    binder_records = [r for r in records if r.chain_id == binder]
    passes = (
        asp6.klass in GATE_ASP6_CLASSES
        and peptide_share >= GATE_MIN_PEPTIDE_SHARE
    )

    sse = backbones.secondary_structure_fraction(array, binder)
    carboxylate = structures.select_atoms(
        array, PEPTIDE_CHAIN, MUTANT_POSITION, ASP_CARBOXYLATE_ATOMS
    )
    binder_coords = array.coord[array.chain_id == binder]
    target_coords = array.coord[np.isin(array.chain_id, TARGET_CHAINS)]

    row = {
        "backbone_id": backbone_id,
        "batch": log_row.get("batch", ""),
        "seed": log_row.get("seed", ""),
        "binder_chain": binder,
        "binder_length": backbones.chain_length(array, binder),
        "target_rmsd": round(
            backbones.target_backbone_rmsd(array, reference, TARGET_CHAINS), 4
        ),
        "asp6_class": asp6.klass,
        "asp6_delta_rsasa": asp6.delta_rsasa,
        "asp6_sasa_loss": round(asp6.sasa_loss, 2),
        "peptide_interface_area": round(peptide_area, 2),
        "hla_interface_area": round(hla_area, 2),
        "peptide_share": round(peptide_share, 4),
        "n_core_binder": sum(1 for r in binder_records if r.klass == "core"),
        "n_support_binder": sum(1 for r in binder_records if r.klass == "support"),
        "n_rim_binder": sum(1 for r in binder_records if r.klass == "rim"),
        "n_interface_binder": sum(
            1 for r in binder_records if r.klass in interface.INTERFACE_CLASSES
        ),
        "n_interface_target": sum(
            1 for r in interface.interface_residues(records) if r.chain_id != binder
        ),
        "buried_interface_area": round(
            sum(r.sasa_loss for r in interface.buried_interface_residues(records)), 2
        ),
        "helix": sse["helix"],
        "strand": sse["strand"],
        "coil": sse["coil"],
        "min_dist_to_target": round(min_distance(binder_coords, target_coords), 2),
        "min_dist_to_carboxylate": round(min_distance(binder_coords, carboxylate), 2),
        "passes_gate": passes,
    }
    return row, records


def write_interface_table(backbone_id: str, records: list) -> None:
    INTERFACE_DIR.mkdir(parents=True, exist_ok=True)
    fields = (
        "chain_id", "res_id", "res_name", "group", "sasa_monomer", "sasa_complex",
        "rsasa_monomer", "rsasa_complex", "delta_rsasa", "changed", "klass",
    )
    with (INTERFACE_DIR / f"{backbone_id}.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow({f: getattr(record, f) for f in fields})


def validate(rows: list[dict]) -> list[str]:
    """Hard checks. These are bugs or wrong inputs, not design failures."""
    problems = []
    for row in rows:
        name = row["backbone_id"]
        if row["target_rmsd"] >= MAX_TARGET_RMSD:
            problems.append(
                f"{name}: target backbone moved {row['target_rmsd']} A - the "
                "output was misparsed or the contig was wrong"
            )
        if not MIN_BINDER_LENGTH <= row["binder_length"] <= MAX_BINDER_LENGTH:
            problems.append(
                f"{name}: binder length {row['binder_length']} outside "
                f"{MIN_BINDER_LENGTH}-{MAX_BINDER_LENGTH}"
            )
        if row["min_dist_to_target"] < MIN_HEAVY_ATOM_APPROACH:
            problems.append(
                f"{name}: closest approach to the target "
                f"{row['min_dist_to_target']} A - the binder overlaps it"
            )
        # Occlusion is mutual, so an interface cannot exist on one side only.
        # When it appears to, the generated chain is inside the target rather
        # than against it, and the per-chain areas are meaningless.
        one_sided = bool(row["n_interface_binder"]) != bool(row["n_interface_target"])
        if one_sided:
            problems.append(
                f"{name}: interface is one-sided "
                f"(binder {row['n_interface_binder']}, "
                f"target {row['n_interface_target']} residues) - the chain is "
                "interpenetrating, not docked"
            )
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--drive-dir", required=True, type=Path)
    parser.add_argument(
        "--target",
        default="7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb",
        help="reference target in data/targets_relaxed/, for the RMSD check",
    )
    parser.add_argument(
        "--relabel-contigs",
        help=(
            "contig string, e.g. 'A1-180/0 C1-10/0 60-75'. Use when the inputs "
            "are raw RFdiffusion output, which is one continuously renumbered "
            "chain: the target's chain ids and numbering are restored here "
            "rather than in the notebook, so there is one tested implementation "
            "of it (kras_pmhc.nim.relabel_chains) instead of one per transport."
        ),
    )
    args = parser.parse_args()

    drive_dir = args.drive_dir.expanduser()
    source_dir = drive_dir / "backbones"
    if not source_dir.is_dir():
        raise SystemExit(f"{source_dir} not found")

    generated = sorted(source_dir.glob("*.pdb"))
    if not generated:
        raise SystemExit(f"no .pdb files in {source_dir}")

    reference = load_reference(args.target)
    log = read_run_log(drive_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"ingesting {len(generated)} backbone(s) from {source_dir}")
    if args.relabel_contigs:
        print(f"relabelling chains from contigs: {args.relabel_contigs}")
    rows = []
    for path in generated:
        backbone_id = path.stem
        local = OUT_DIR / path.name
        shutil.copy2(path, local)

        array = read_pdb(local)
        if args.relabel_contigs:
            # Raw RFdiffusion output is one continuously renumbered chain. The
            # corrected file is written back, so every later stage reads chains
            # A/C/B and no stage has to repeat this.
            array = nim.relabel_chains(array, args.relabel_contigs)
            relabelled = PDBFile()
            relabelled.set_structure(array)
            relabelled.write(str(local))

        row, records = score(array, reference, backbone_id, log.get(path.name, {}))
        write_interface_table(backbone_id, records)
        rows.append(row)

        flag = "PASS" if row["passes_gate"] else "fail"
        print(
            f"  {backbone_id:<16} {row['binder_length']:>3} aa  "
            f"Asp6 {row['asp6_class']:<8} peptide share "
            f"{row['peptide_share']:.0%}  {flag}"
        )

    manifest = OUT_DIR / "MANIFEST.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    problems = validate(rows)
    passing = [r for r in rows if r["passes_gate"]]

    print(f"\nwrote {manifest.relative_to(ROOT)}")
    print(f"      {INTERFACE_DIR.relative_to(ROOT)}/<backbone_id>.tsv")
    print(f"\npilot gate: {len(passing)}/{len(rows)} engage Asp6 "
          f"(need >= {GATE_MIN_PASSING} of 8 to proceed)")
    print(f"  Asp6 is support/rim/core and peptide carries "
          f">= {GATE_MIN_PEPTIDE_SHARE:.0%} of target interface area")

    if problems:
        print("\nVALIDATION FAILED:")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    if len(rows) >= GATE_MIN_PASSING and not passing:
        print(
            "\nNo backbone engages the peptide. The binder is landing on the "
            "alpha1/alpha2 helix tops.\nRemedies, one pilot each: widen "
            "hotspots to [C4,C5,C6,C7,C8]; then use_beta_model."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
