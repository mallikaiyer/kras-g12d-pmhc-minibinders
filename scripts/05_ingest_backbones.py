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
#: Shape of the target after superposition, over N/CA/C. This is a parsing
#: check rather than physics, and the two outcomes it separates are far apart:
#: the target's frames are reproduced to ~0.1 A, while a mis-assigned chain
#: ordering gave 21.6 A. 0.5 sits between them with room for both.
MAX_TARGET_RMSD = 0.5
MIN_BINDER_LENGTH = 60
MAX_BINDER_LENGTH = 75
#: Below this the generated chain overlaps the target rather than touching it.
#: Measured against the *whole* target, not the carboxylate: a binder can sit a
#: comfortable distance from Asp6 while interpenetrating the HLA elsewhere, and
#: that failure is invisible to a carboxylate-only check.
MIN_HEAVY_ATOM_APPROACH = 2.0

MANIFEST_FIELDS = (
    "backbone_id", "batch", "seed", "binder_chain", "binder_length",
    "target_rmsd", "target_offset",
    "asp6_class", "asp6_delta_rsasa", "asp6_sasa_loss",
    "peptide_interface_area", "hla_interface_area", "peptide_share",
    "n_core_binder", "n_support_binder", "n_rim_binder",
    "n_interface_binder", "n_interface_target", "buried_interface_area",
    "helix", "strand", "coil",
    "min_dist_to_target", "min_dist_to_carboxylate", "clashes", "passes_gate",
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


def score(
    array, shape_rmsd: float, offset: float, backbone_id: str, log_row: dict
) -> tuple[dict, list]:
    """Metrics for one generated backbone, plus its per-residue Levy table.

    ``array`` is the spliced complex - the full-atom stage-04 target with the
    generated binder added - not raw generator output. The interface numbers
    would be meaningless otherwise: RFdiffusion emits backbone atoms only, and
    a target without side chains has nothing to bury.
    """
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

    sse = backbones.secondary_structure_fraction(array, binder)
    carboxylate = structures.select_atoms(
        array, PEPTIDE_CHAIN, MUTANT_POSITION, ASP_CARBOXYLATE_ATOMS
    )
    binder_coords = array.coord[array.chain_id == binder]
    target_coords = array.coord[np.isin(array.chain_id, TARGET_CHAINS)]

    # Measured against the whole target, side chains included: a binder can sit
    # a comfortable distance from Asp6 while overlapping the HLA elsewhere.
    closest = min_distance(binder_coords, target_coords)
    clashes = closest < MIN_HEAVY_ATOM_APPROACH

    # A clashing design cannot engage anything - its interface areas are
    # computed from interpenetrating atoms - so it fails the gate whatever the
    # Levy classification says.
    passes = (
        not clashes
        and asp6.klass in GATE_ASP6_CLASSES
        and peptide_share >= GATE_MIN_PEPTIDE_SHARE
    )

    row = {
        "backbone_id": backbone_id,
        "batch": log_row.get("batch", ""),
        "seed": log_row.get("seed", ""),
        "binder_chain": binder,
        "binder_length": backbones.chain_length(array, binder),
        # Shape check, after superposition on the target backbone.
        "target_rmsd": round(shape_rmsd, 4),
        # Rigid-body displacement of the generator's output from the reference.
        # Reported, not judged: RFdiffusion recentres, so tens of angstroms
        # here is normal and says nothing about whether the geometry survived.
        "target_offset": round(offset, 2),
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
        "min_dist_to_target": round(closest, 2),
        "clashes": clashes,
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
    """Hard checks: bugs or wrong inputs, not design failures.

    A steric clash is deliberately *not* here. RFdiffusion sometimes places a
    backbone inside its target, and that is a property of the design, not of
    this pipeline - at 50 designs several did. Clashing designs are recorded in
    the manifest, excluded from the gate, and reported, but they do not stop
    the run; aborting on one would discard the other 49.
    """
    problems = []
    for row in rows:
        name = row["backbone_id"]
        if row["target_rmsd"] >= MAX_TARGET_RMSD:
            problems.append(
                f"{name}: target shape changed by {row['target_rmsd']} A "
                "after superposition - the output was misparsed, the contig "
                "was wrong, or the chains were assigned in the wrong order"
            )
        if not MIN_BINDER_LENGTH <= row["binder_length"] <= MAX_BINDER_LENGTH:
            problems.append(
                f"{name}: binder length {row['binder_length']} outside "
                f"{MIN_BINDER_LENGTH}-{MAX_BINDER_LENGTH}"
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
            # Raw RFdiffusion output is one continuously renumbered chain.
            array = nim.relabel_chains(array, args.relabel_contigs)

        # Generator output is backbone-only and recentred. Undo the recentering
        # by superimposing on the target, then reattach the target's real side
        # chains by splicing the binder onto the stage-04 reference. The result
        # is what gets stored and scored: every later stage wants the full-atom
        # complex in the stage-04 frame, not raw generator output.
        binder = backbones.binder_chain(array, TARGET_CHAINS)
        offset = backbones.target_backbone_rmsd(array, reference, TARGET_CHAINS)
        aligned, shape_rmsd = backbones.superimpose_on_target(
            array, reference, TARGET_CHAINS
        )
        complexed = backbones.splice_binder(reference, aligned, binder)

        spliced = PDBFile()
        spliced.set_structure(complexed)
        spliced.write(str(local))

        row, records = score(
            complexed, shape_rmsd, offset, backbone_id, log.get(path.name, {})
        )
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
    clashing = [row for row in rows if row["clashes"]]

    print(f"\ngate: {len(passing)}/{len(rows)} engage Asp6")
    print(f"  Asp6 is support/rim/core and peptide carries "
          f">= {GATE_MIN_PEPTIDE_SHARE:.0%} of target interface area")
    if clashing:
        print(f"  {len(clashing)} rejected for overlapping the target "
              f"(< {MIN_HEAVY_ATOM_APPROACH} A): "
              + ", ".join(row["backbone_id"] for row in clashing))
    if len(rows) == 8 and len(passing) < GATE_MIN_PASSING:
        print(f"  a pilot of 8 needs >= {GATE_MIN_PASSING} passing to justify "
              "a production run")

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
