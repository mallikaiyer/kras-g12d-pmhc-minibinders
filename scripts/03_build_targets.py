#!/usr/bin/env python
"""Stage 03 - build the Route A design targets.

Takes the two complete wild-type copies from 7OW3, whose peptide position 6
projects out of the HLA groove, and mutates Gly6 to Asp at the rotamers
identified by the stage-02 scan. The backbone is left exactly as deposited:
only four side-chain atoms are added.

For each (copy, rotamer) this writes:

  <id>.pdb        full pMHC - peptide, heavy chain, beta-2-microglobulin.
                  Input to stage 04 (relaxation and MD) and the co-folding
                  reference.
  <id>_crop.pdb   HLA alpha1+alpha2 plus peptide only. Input to stage 05
                  generation; beta-2-microglobulin and alpha3 sit on the far
                  side of the platform and only cost GPU memory.

Plus, once per copy:

  <copy>_WT.pdb   the unmutated complex, for the stage-07 differential.
  <copy>_apo.pdb  peptide deleted, the objective-2 counter-target.

Chains are renamed to the pMHC convention A=heavy, B=beta-2-microglobulin,
C=peptide so that stage 05 contigs are identical across targets. The original
chain letters are recorded in the manifest.

Every written file is read back and checked before the script exits.

Usage:
    python scripts/03_build_targets.py
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

from kras_pmhc import structures, survey  # noqa: E402
from kras_pmhc.build import mutate_gly_to_asp  # noqa: E402
from kras_pmhc.constants import (  # noqa: E402
    GROOVE_REFERENCE_ATOM,
    MUTANT_POSITION,
    PEPTIDE_G12D,
)
from kras_pmhc.geometry import min_distance  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "targets"
SOURCE_PDB = "7OW3"

#: Canonical chain naming for the written targets.
CHAIN_HEAVY, CHAIN_B2M, CHAIN_PEPTIDE = "A", "B", "C"

#: Last heavy-chain residue kept in the cropped generation target. The class I
#: heavy chain is alpha1 (1-90), alpha2 (91-182), alpha3 (183-275); 1-180
#: captures the peptide-binding platform and its two flanking helices.
CROP_LAST_RESIDUE = 180

#: Rotamers to build, chosen in stage 02 as the best-exposed clash-free
#: rotamer in each canonical chi1 well. Copy L has no viable +60 rotamer.
#: (peptide_chain, chi1, chi2)
ROTAMERS = [
    ("C", -70.0, 120.0),
    ("C", -180.0, 30.0),
    ("C", 60.0, 45.0),
    ("L", -70.0, 120.0),
    ("L", -180.0, 0.0),
]

#: Exposure a built carboxylate must retain to be worth designing against.
#: The closed state measures ~7 A^2.
MIN_CARBOXYLATE_SASA = 25.0

CARBOXYLATE = ("CG", "OD1", "OD2")
SIDECHAIN = ("CB", "CG", "OD1", "OD2")


def angle_tag(value: float) -> str:
    """Filesystem-safe label for a torsion, e.g. -70 -> 'm70', 120 -> 'p120'."""
    return f"{'m' if value < 0 else 'p'}{abs(int(round(value)))}"


def rename_chains(array: struc.AtomArray, peptide: str, heavy: str, b2m: str):
    """Relabel to the A/B/C pMHC convention, leaving coordinates untouched."""
    out = array.copy()
    original = out.chain_id.copy()
    out.chain_id[original == heavy] = CHAIN_HEAVY
    out.chain_id[original == b2m] = CHAIN_B2M
    out.chain_id[original == peptide] = CHAIN_PEPTIDE
    return out


def write_pdb(array: struc.AtomArray, path: Path) -> None:
    pdb = PDBFile()
    pdb.set_structure(array)
    pdb.write(str(path))


def torsions_at(array: struc.AtomArray, chain: str, res_id: int) -> tuple[float, float]:
    """Measure chi1 and chi2 of an aspartate, in degrees."""
    get = lambda name: structures.atom_coord(array, chain, res_id, name)  # noqa: E731
    chi1 = float(np.rad2deg(struc.dihedral(get("N"), get("CA"), get("CB"), get("CG"))))
    chi2 = float(np.rad2deg(struc.dihedral(get("CA"), get("CB"), get("CG"), get("OD1"))))
    return chi1, chi2


def subset_sasa(array: struc.AtomArray, chain: str, res_id: int, names) -> float:
    """SASA summed over named atoms of one residue, in the context of ``array``."""
    mask = (
        (array.chain_id == chain)
        & (array.res_id == res_id)
        & np.isin(array.atom_name, list(names))
    )
    if not mask.any():
        raise ValueError(f"no atoms {tuple(names)} in chain {chain} residue {res_id}")
    return float(np.nansum(struc.sasa(array)[mask]))


def report_gaps(array: struc.AtomArray, chain: str, last: int) -> str:
    """List residue numbers missing from 1..last in a chain."""
    present = set(int(r) for r in np.unique(array.res_id[array.chain_id == chain]))
    missing = [r for r in range(1, last + 1) if r not in present]
    return ",".join(str(m) for m in missing) if missing else ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    source = structures.load(SOURCE_PDB)
    peptide_to_heavy = structures.pair_peptides_with_heavy_chains(source)
    heavy_to_b2m = structures.pair_b2m_with_heavy_chains(source)
    ref_resname, ref_res_id, ref_atom = GROOVE_REFERENCE_ATOM

    # Reference phi/psi of the untouched glycine, to prove the mutation does
    # not move the backbone.
    wt_dihedrals = {
        peptide: survey.phi_psi(source, peptide, MUTANT_POSITION)
        for peptide, _, _ in ROTAMERS
    }

    rows: list[dict] = []
    copies_done: set[str] = set()

    for peptide, chi1, chi2 in ROTAMERS:
        heavy = peptide_to_heavy[peptide]
        b2m = heavy_to_b2m[heavy]
        complex_wt = structures.extract_complex(source, peptide)

        # --- the unmutated and apo references, once per copy ---------------
        if peptide not in copies_done:
            copies_done.add(peptide)
            renamed_wt = rename_chains(complex_wt, peptide, heavy, b2m)
            write_pdb(renamed_wt, args.out_dir / f"7ow3_{peptide}_WT.pdb")
            write_pdb(
                renamed_wt[renamed_wt.chain_id != CHAIN_PEPTIDE],
                args.out_dir / f"7ow3_{peptide}_apo.pdb",
            )

        # --- the mutant ----------------------------------------------------
        mutated = mutate_gly_to_asp(complex_wt, peptide, MUTANT_POSITION, chi1, chi2)
        renamed = rename_chains(mutated, peptide, heavy, b2m)

        target_id = f"7ow3_{peptide}_G6D_chi1{angle_tag(chi1)}_chi2{angle_tag(chi2)}"
        full_path = args.out_dir / f"{target_id}.pdb"
        crop_path = args.out_dir / f"{target_id}_crop.pdb"
        write_pdb(renamed, full_path)

        crop = renamed[
            ((renamed.chain_id == CHAIN_HEAVY) & (renamed.res_id <= CROP_LAST_RESIDUE))
            | (renamed.chain_id == CHAIN_PEPTIDE)
        ]
        write_pdb(crop, crop_path)

        # --- read back and validate ---------------------------------------
        check = PDBFile.read(str(full_path)).get_structure(
            model=1, extra_fields=["b_factor", "occupancy"]
        )

        sequence = structures.observed_sequence(check, CHAIN_PEPTIDE, len(PEPTIDE_G12D))
        if sequence != PEPTIDE_G12D:
            raise ValueError(f"{target_id}: peptide is {sequence}, expected {PEPTIDE_G12D}")

        n_chains = len(np.unique(check.chain_id))
        if n_chains != 3:
            raise ValueError(f"{target_id}: {n_chains} chains written, expected 3")

        got_chi1, got_chi2 = torsions_at(check, CHAIN_PEPTIDE, MUTANT_POSITION)
        for label, want, got in (("chi1", chi1, got_chi1), ("chi2", chi2, got_chi2)):
            # Torsions are circular; compare on the shortest arc.
            delta = abs((want - got + 180.0) % 360.0 - 180.0)
            if delta > 0.5:
                raise ValueError(f"{target_id}: {label} is {got:.1f}, requested {want:.1f}")

        phi, psi = survey.phi_psi(check, CHAIN_PEPTIDE, MUTANT_POSITION)
        want_phi, want_psi = wt_dihedrals[peptide]
        for label, want, got in (("phi", want_phi, phi), ("psi", want_psi, psi)):
            if abs(want - got) > 0.1:
                raise ValueError(
                    f"{target_id}: {label} moved from {want:.1f} to {got:.1f}; "
                    "the mutation must not perturb the backbone"
                )

        ref = structures.atom_coord(check, CHAIN_HEAVY, ref_res_id, ref_atom)
        carboxylate = structures.select_atoms(
            check, CHAIN_PEPTIDE, MUTANT_POSITION, CARBOXYLATE
        )
        carbox_to_ref = min_distance(carboxylate, np.array([ref]))

        others = check[
            ~((check.chain_id == CHAIN_PEPTIDE) & (check.res_id == MUTANT_POSITION))
        ]
        closest = min_distance(
            structures.select_atoms(check, CHAIN_PEPTIDE, MUTANT_POSITION, SIDECHAIN),
            others.coord,
        )

        carbox_sasa = subset_sasa(check, CHAIN_PEPTIDE, MUTANT_POSITION, CARBOXYLATE)
        sidechain_sasa = subset_sasa(check, CHAIN_PEPTIDE, MUTANT_POSITION, SIDECHAIN)
        if carbox_sasa < MIN_CARBOXYLATE_SASA:
            raise ValueError(
                f"{target_id}: carboxylate SASA {carbox_sasa:.1f} A^2 is below the "
                f"{MIN_CARBOXYLATE_SASA} A^2 floor; this rotamer is not engageable"
            )

        rows.append(
            {
                "target_id": target_id,
                "source": SOURCE_PDB,
                "source_peptide_chain": peptide,
                "source_heavy_chain": heavy,
                "source_b2m_chain": b2m,
                "chi1": round(chi1, 1),
                "chi2": round(chi2, 1),
                "phi": round(phi, 1),
                "psi": round(psi, 1),
                "carboxylate_sasa": round(carbox_sasa, 2),
                "sidechain_sasa": round(sidechain_sasa, 2),
                "carboxylate_to_ref": round(carbox_to_ref, 2),
                "closest_heavy_atom": round(closest, 2),
                "reference": f"{ref_resname}{ref_res_id}:{ref_atom}",
                "n_atoms_full": len(check),
                "n_atoms_crop": len(crop),
                "heavy_gaps_in_crop": report_gaps(crop, CHAIN_HEAVY, CROP_LAST_RESIDUE),
                "full_pdb": full_path.name,
                "crop_pdb": crop_path.name,
            }
        )
        print(
            f"  {target_id}\n"
            f"      chi1={got_chi1:7.1f} chi2={got_chi2:7.1f}  "
            f"carboxylate SASA {carbox_sasa:5.1f} A^2  "
            f"-> {ref_resname}{ref_res_id} {carbox_to_ref:5.2f} A  "
            f"closest {closest:4.2f} A"
        )

    manifest = args.out_dir / "MANIFEST.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n{len(rows)} mutant targets (+ WT and apo per copy) in "
          f"{args.out_dir.relative_to(REPO_ROOT)}")
    print(f"manifest: {manifest.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
