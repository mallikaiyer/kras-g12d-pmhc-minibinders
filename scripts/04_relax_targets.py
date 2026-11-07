#!/usr/bin/env python
"""Stage 04 - protonate and relax the design targets.

Adds hydrogens and runs a restrained energy minimisation on every structure
stage 03 produced: five mutants, two wild-type references and two apo
references. All nine get an identical protocol, which matters because the
wild-type reference is the comparator for the stage-07 differential - if it
were left unrelaxed it would differ from the mutant for reasons unrelated to
the mutation.

Backbone heavy atoms are restrained; every side chain is free. This preserves
the experimentally observed backbone while letting the environment respond to
the new aspartate.

**This stage runs no dynamics.** It says nothing about whether the open
conformer is populated in solution. That is a stated shortcoming of this
submission, not an oversight - see docs/design_plan.md.

Usage:
    python scripts/04_relax_targets.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import biotite.structure as struc
import numpy as np
from biotite.structure.io.pdb import PDBFile

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from kras_pmhc import relax, structures, survey  # noqa: E402
from kras_pmhc.constants import (  # noqa: E402
    GROOVE_REFERENCE_ATOM,
    HLA_SURVEY_RESIDUES,
    MUTANT_POSITION,
    PEPTIDE_G12D,
)
from kras_pmhc.geometry import min_distance  # noqa: E402

IN_DIR = REPO_ROOT / "data" / "targets"
OUT_DIR = REPO_ROOT / "data" / "targets_relaxed"

CHAIN_HEAVY, CHAIN_B2M, CHAIN_PEPTIDE = "A", "B", "C"
CROP_LAST_RESIDUE = 180
BACKBONE = ("N", "CA", "C", "O")
CARBOXYLATE = ("CG", "OD1", "OD2")
ASP_SIDECHAIN = ("CB", "CG", "OD1", "OD2")

#: Groove side chains whose response to the new aspartate is worth reporting.
#: If minimisation drags Arg114 up toward an open Asp6, that is evidence of
#: residual pull toward the closed state.
WATCHED_HLA_RESIDUES = (70, 114, 155)

# --- validation thresholds -----------------------------------------------
MAX_BACKBONE_RMSD = 0.3        # A; the restraint should hold the trace
# Dihedrals amplify small positional changes: a 0.1 A backbone RMSD routinely
# shifts phi/psi by a few degrees. The backbone RMSD bound above is the real
# guard on the trace, so this only needs to catch a genuine excursion.
MAX_PHI_PSI_DRIFT = 15.0       # degrees
# chi1 places the carboxylate in space, so a large change means the rotamer
# well was not a minimum. chi2 merely spins a near-symmetric planar carboxylate
# about the CB-CG axis; it is a shallow, broadly distributed torsion and tens of
# degrees of movement carries no information. Gate the first, report the second.
MAX_CHI1_DRIFT = 30.0          # degrees
MAX_CHI2_DRIFT = 90.0          # degrees; a sanity bound only
MIN_CARBOXYLATE_SASA = 25.0    # A^2
MIN_HEAVY_ATOM_APPROACH = 2.4  # A


def read_heavy(path: Path) -> struc.AtomArray:
    """Read a PDB and drop hydrogens."""
    array = PDBFile.read(str(path)).get_structure(model=1)
    return array[array.element != "H"]


def coord_map(array: struc.AtomArray) -> dict[tuple[str, int, str], np.ndarray]:
    return {
        (str(array.chain_id[i]), int(array.res_id[i]), str(array.atom_name[i])): array.coord[i]
        for i in range(len(array))
    }


def rmsd_over(
    before: dict, after: dict, predicate
) -> tuple[float | None, int]:
    """RMSD over atoms selected by ``predicate(chain, res_id, atom_name)``.

    Only atoms present in both structures are compared. Returns ``(None, 0)``
    when nothing matches, which callers treat as "not applicable" rather than
    as zero displacement.
    """
    deltas = [
        float(np.linalg.norm(after[key] - coord))
        for key, coord in before.items()
        if key in after and predicate(*key)
    ]
    if not deltas:
        return None, 0
    return float(np.sqrt(np.mean(np.square(deltas)))), len(deltas)


def torsions(array: struc.AtomArray) -> tuple[float, float] | tuple[None, None]:
    get = lambda n: structures.atom_coord(array, CHAIN_PEPTIDE, MUTANT_POSITION, n)  # noqa: E731
    try:
        n, ca, cb, cg, od1 = (get(x) for x in ("N", "CA", "CB", "CG", "OD1"))
    except Exception:
        return (None, None)
    if any(x is None for x in (n, ca, cb, cg, od1)):
        return (None, None)
    return (
        float(np.rad2deg(struc.dihedral(n, ca, cb, cg))),
        float(np.rad2deg(struc.dihedral(ca, cb, cg, od1))),
    )


def circular_delta(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def subset_sasa(array: struc.AtomArray, names) -> float | None:
    mask = (
        (array.chain_id == CHAIN_PEPTIDE)
        & (array.res_id == MUTANT_POSITION)
        & np.isin(array.atom_name, list(names))
    )
    if not mask.any():
        return None
    return float(np.nansum(struc.sasa(array)[mask]))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-dir", type=Path, default=IN_DIR)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--restraint-k", type=float, default=relax.DEFAULT_RESTRAINT_K)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sources = sorted(p for p in args.in_dir.glob("*.pdb") if not p.stem.endswith("_crop"))
    if not sources:
        print(f"no targets in {args.in_dir}; run stage 03 first", file=sys.stderr)
        return 1

    ref_resname, ref_res_id, ref_atom = GROOVE_REFERENCE_ATOM
    rows: list[dict] = []

    for source in sources:
        target_id = source.stem
        before = read_heavy(source)
        is_mutant = "G6D" in target_id
        has_peptide = CHAIN_PEPTIDE in set(map(str, before.chain_id))

        topology, positions = relax.protonate(source)
        relaxed, energy_before, energy_after, n_restrained = relax.minimise(
            topology, positions, restraint_k=args.restraint_k
        )

        full = args.out_dir / f"{target_id}_relaxed.pdb"
        relax.write_pdb(topology, relaxed, full)

        after = read_heavy(full)
        heavy_path = args.out_dir / f"{target_id}_relaxed_heavy.pdb"
        pdb = PDBFile()
        pdb.set_structure(after)
        pdb.write(str(heavy_path))

        crop_path = None
        if has_peptide:
            crop = after[
                ((after.chain_id == CHAIN_HEAVY) & (after.res_id <= CROP_LAST_RESIDUE))
                | (after.chain_id == CHAIN_PEPTIDE)
            ]
            crop_path = args.out_dir / f"{target_id}_relaxed_crop.pdb"
            pdb = PDBFile()
            pdb.set_structure(crop)
            pdb.write(str(crop_path))

        # --- metrics -----------------------------------------------------
        b_map, a_map = coord_map(before), coord_map(after)
        backbone_rmsd, n_bb = rmsd_over(b_map, a_map, lambda c, r, n: n in BACKBONE)
        all_heavy_rmsd, _ = rmsd_over(b_map, a_map, lambda c, r, n: True)
        asp_rmsd, _ = rmsd_over(
            b_map,
            a_map,
            lambda c, r, n: c == CHAIN_PEPTIDE and r == MUTANT_POSITION and n in ASP_SIDECHAIN,
        )
        watched = {}
        for res_id in WATCHED_HLA_RESIDUES:
            value, _ = rmsd_over(
                b_map,
                a_map,
                lambda c, r, n, _r=res_id: (
                    c == CHAIN_HEAVY and r == _r and n not in BACKBONE
                ),
            )
            label = f"{HLA_SURVEY_RESIDUES.get(res_id, '')}{res_id}"
            watched[label] = round(value, 3) if value is not None else None

        chi1_before, chi2_before = torsions(before)
        chi1_after, chi2_after = torsions(after)
        phi_before, psi_before = (
            survey.phi_psi(before, CHAIN_PEPTIDE, MUTANT_POSITION) if has_peptide else (None, None)
        )
        phi_after, psi_after = (
            survey.phi_psi(after, CHAIN_PEPTIDE, MUTANT_POSITION) if has_peptide else (None, None)
        )

        carbox_sasa = subset_sasa(after, CARBOXYLATE)
        carbox_to_ref = None
        closest = None
        if is_mutant:
            ref = structures.atom_coord(after, CHAIN_HEAVY, ref_res_id, ref_atom)
            carbox_to_ref = min_distance(
                structures.select_atoms(after, CHAIN_PEPTIDE, MUTANT_POSITION, CARBOXYLATE),
                np.array([ref]),
            )
            others = after[
                ~((after.chain_id == CHAIN_PEPTIDE) & (after.res_id == MUTANT_POSITION))
            ]
            closest = min_distance(
                structures.select_atoms(after, CHAIN_PEPTIDE, MUTANT_POSITION, ASP_SIDECHAIN),
                others.coord,
            )

        # --- validation --------------------------------------------------
        if backbone_rmsd is None or backbone_rmsd > MAX_BACKBONE_RMSD:
            raise ValueError(
                f"{target_id}: backbone RMSD {backbone_rmsd} A exceeds "
                f"{MAX_BACKBONE_RMSD} A; the restraint is not holding"
            )
        if has_peptide:
            sequence = structures.observed_sequence(after, CHAIN_PEPTIDE, len(PEPTIDE_G12D))
            expected = PEPTIDE_G12D if is_mutant else "VVVGAGGVGK"
            if sequence != expected:
                raise ValueError(f"{target_id}: peptide {sequence}, expected {expected}")
            for label, b, a in (("phi", phi_before, phi_after), ("psi", psi_before, psi_after)):
                if b is not None and a is not None and circular_delta(b, a) > MAX_PHI_PSI_DRIFT:
                    raise ValueError(
                        f"{target_id}: {label} drifted {circular_delta(b, a):.1f} deg "
                        f"({b:.1f} -> {a:.1f})"
                    )
            # The substantive Ramachandran claim from stage 02: position 6 has
            # negative phi, so the backbone is not glycine-exclusive and an
            # aspartate can occupy it. Relaxation must not overturn that.
            if phi_after is not None and phi_after >= 0.0:
                raise ValueError(
                    f"{target_id}: phi relaxed to {phi_after:.1f} deg (positive); "
                    "the open backbone would no longer admit a non-glycine residue"
                )
        if is_mutant:
            for label, b, a, limit in (
                ("chi1", chi1_before, chi1_after, MAX_CHI1_DRIFT),
                ("chi2", chi2_before, chi2_after, MAX_CHI2_DRIFT),
            ):
                if circular_delta(b, a) > limit:
                    raise ValueError(
                        f"{target_id}: {label} moved {circular_delta(b, a):.1f} deg "
                        f"({b:.1f} -> {a:.1f}), past its {limit:.0f} deg bound"
                    )
            if carbox_sasa < MIN_CARBOXYLATE_SASA:
                raise ValueError(
                    f"{target_id}: carboxylate SASA {carbox_sasa:.1f} A^2 below floor"
                )
            if closest < MIN_HEAVY_ATOM_APPROACH:
                raise ValueError(f"{target_id}: clash at {closest:.2f} A")

        rows.append(
            {
                "target_id": target_id,
                "is_mutant": is_mutant,
                "n_atoms_heavy": len(after),
                "n_restrained_backbone": n_restrained,
                "energy_before_kcal": round(energy_before, 1),
                "energy_after_kcal": round(energy_after, 1),
                "backbone_rmsd": round(backbone_rmsd, 3),
                "all_heavy_rmsd": round(all_heavy_rmsd, 3),
                "asp6_sidechain_rmsd": round(asp_rmsd, 3) if asp_rmsd is not None else None,
                "chi1_before": round(chi1_before, 1) if chi1_before is not None else None,
                "chi1_after": round(chi1_after, 1) if chi1_after is not None else None,
                "chi2_before": round(chi2_before, 1) if chi2_before is not None else None,
                "chi2_after": round(chi2_after, 1) if chi2_after is not None else None,
                "phi_after": round(phi_after, 1) if phi_after is not None else None,
                "psi_after": round(psi_after, 1) if psi_after is not None else None,
                "carboxylate_sasa": round(carbox_sasa, 2) if carbox_sasa is not None else None,
                "carboxylate_to_ref": (
                    round(carbox_to_ref, 2) if carbox_to_ref is not None else None
                ),
                "closest_heavy_atom": round(closest, 2) if closest is not None else None,
                **{f"sidechain_rmsd_{k}": v for k, v in watched.items()},
                "relaxed_pdb": full.name,
                "heavy_pdb": heavy_path.name,
                "crop_pdb": crop_path.name if crop_path else "",
            }
        )

        asp = f"chi1 {chi1_before:7.1f}->{chi1_after:7.1f}" if is_mutant else "  (no Asp6)       "
        print(
            f"  {target_id:36s} bbRMSD {backbone_rmsd:.3f}  {asp}  "
            f"E {energy_before:10.1f} -> {energy_after:9.1f} kcal/mol"
        )

    manifest = args.out_dir / "MANIFEST.tsv"
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if row.get(k) is None else row.get(k) for k in fieldnames})

    print(f"\n{len(rows)} structures relaxed into {args.out_dir.relative_to(REPO_ROOT)}")
    print(f"manifest: {manifest.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
