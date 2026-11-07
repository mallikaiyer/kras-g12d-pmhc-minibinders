#!/usr/bin/env python
"""Stage 02 - can an aspartate on the open backbone actually be engaged?

The position-6 orientation table in the survey measures where Cbeta sits.
Cbeta is fully determined by the backbone, but the carboxylate is not: chi1 and
chi2 are free, so that table bounds the carboxylate's position without fixing
its direction. This stage resolves the remaining freedom explicitly.

For each complete wild-type copy in 7OW3, Gly6 is replaced by an aspartate
built at every (chi1, chi2) on a grid. Each rotamer is then scored for:

  * steric viability - closest approach to any other heavy atom in the complex
  * reachability     - distance from the carboxylate to HLA Arg114 CZ, the
                       partner it pairs with in the closed state
  * exposure         - solvent-accessible area of the built side chain

The question this answers: of the rotamers an Asp could adopt on this
backbone, how many are clash-free, and are those exposed enough for a designed
binder to engage?

Usage:
    python scripts/02_asp_rotamer_scan.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import biotite.structure as struc
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from kras_pmhc import structures  # noqa: E402
from kras_pmhc.build import ASP_CHI1_WELLS, ASP_CHI2_PERIOD, build_asp_sidechain  # noqa: E402
from kras_pmhc.constants import GROOVE_REFERENCE_ATOM, MUTANT_POSITION  # noqa: E402
from kras_pmhc.geometry import min_distance, virtual_cb  # noqa: E402

OUT_DIR = REPO_ROOT / "results" / "rotamer_scan"

#: Any built heavy atom closer than this to another heavy atom is a hard
#: clash. Kept below hydrogen-bond range (2.6-3.2 A) so that a genuine
#: carboxylate-to-donor contact is not misreported as a collision.
CLASH_CUTOFF = 2.4

#: Side-chain exposure above which a carboxylate is plausibly engageable by a
#: designed binder. A fully buried Asp in the closed state measures ~7 A^2.
EXPOSURE_THRESHOLD = 25.0


def _append_sidechain(
    array: struc.AtomArray, chain_id: str, res_id: int, atoms: dict[str, np.ndarray]
) -> tuple[struc.AtomArray, np.ndarray]:
    """Return a copy of ``array`` with built atoms appended, and their mask.

    Appending rather than editing in place keeps the scan free of state: each
    rotamer is scored against a pristine structure.
    """
    built = struc.AtomArray(len(atoms))
    for i, (name, coord) in enumerate(atoms.items()):
        built.coord[i] = coord
        built.chain_id[i] = chain_id
        built.res_id[i] = res_id
        built.res_name[i] = "ASP"
        built.atom_name[i] = name
        built.element[i] = name[0]
    combined = array + built
    mask = np.zeros(len(combined), dtype=bool)
    mask[len(array):] = True
    return combined, mask


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb-id", default="7OW3")
    parser.add_argument("--chi1-step", type=float, default=10.0)
    parser.add_argument("--chi2-step", type=float, default=15.0)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    array = structures.load(args.pdb_id)
    pairs = structures.pair_peptides_with_heavy_chains(array)
    ref_resname, ref_res_id, ref_atom = GROOVE_REFERENCE_ATOM

    chi1_grid = np.arange(-180.0, 180.0, args.chi1_step)
    chi2_grid = np.arange(0.0, ASP_CHI2_PERIOD, args.chi2_step)

    rows: list[dict] = []
    for pep, heavy in sorted(pairs.items()):
        res = structures.residue(array, pep, MUTANT_POSITION)
        if len(res) == 0:
            print(f"  {args.pdb_id} chain {pep}: position 6 not modelled, skipping")
            continue
        if str(res.res_name[0]) != "GLY":
            print(f"  {args.pdb_id} chain {pep}: position 6 is {res.res_name[0]}, skipping")
            continue

        n = structures.atom_coord(array, pep, MUTANT_POSITION, "N")
        ca = structures.atom_coord(array, pep, MUTANT_POSITION, "CA")
        c = structures.atom_coord(array, pep, MUTANT_POSITION, "C")
        cb = virtual_cb(n, ca, c)
        ref = structures.atom_coord(array, heavy, ref_res_id, ref_atom)

        # One copy of the complex only: peptide, its heavy chain, and the
        # beta-2-microglobulin paired with it. Neighbouring copies in the
        # asymmetric unit are crystal packing, not biology.
        b2m = [
            c_id
            for c_id, n_res in structures.chain_residue_counts(array).items()
            if 85 <= n_res <= 125
        ]
        context_chains = [pep, heavy] + b2m[:1]
        context = array[np.isin(array.chain_id, context_chains)]
        # Everything except position 6 itself, for clash detection.
        others = context[
            ~((context.chain_id == pep) & (context.res_id == MUTANT_POSITION))
        ]

        print(f"  {args.pdb_id} chain {pep} (heavy {heavy}): scanning "
              f"{len(chi1_grid)}x{len(chi2_grid)} rotamers")

        for chi1 in chi1_grid:
            for chi2 in chi2_grid:
                built = build_asp_sidechain(n, ca, cb, float(chi1), float(chi2))
                coords = np.array(list(built.values()))
                closest = min_distance(coords, others.coord)
                clash = closest < CLASH_CUTOFF

                row = {
                    "pdb_id": args.pdb_id.upper(),
                    "peptide_chain": pep,
                    "chi1": round(float(chi1), 1),
                    "chi2": round(float(chi2), 1),
                    "canonical_chi1_well": min(
                        abs(((chi1 - w + 180) % 360) - 180) for w in ASP_CHI1_WELLS
                    )
                    <= args.chi1_step / 2,
                    "closest_heavy_atom": round(closest, 2),
                    "clash": clash,
                    "carboxylate_to_ref": round(
                        min_distance(
                            np.array([built["OD1"], built["OD2"]]), np.array([ref])
                        ),
                        2,
                    ),
                    "reference": f"{ref_resname}{ref_res_id}:{ref_atom}",
                }

                # SASA is the expensive step; only clash-free rotamers matter.
                if clash:
                    row["sidechain_sasa"] = None
                else:
                    combined, mask = _append_sidechain(
                        context, pep, MUTANT_POSITION, built
                    )
                    per_atom = struc.sasa(combined)
                    # CG + OD1 + OD2 only. The virtual Cbeta is not inserted
                    # into the structure, so it does not occlude its own
                    # carboxylate: exposure here is a slight over-estimate.
                    row["sidechain_sasa"] = round(float(np.nansum(per_atom[mask])), 2)
                rows.append(row)

    out = args.out_dir / "asp_rotamer_scan.tsv"
    fieldnames = list(rows[0])
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if row[k] is None else row[k] for k in fieldnames})

    # --- summary ---------------------------------------------------------
    print()
    for pep in sorted({r["peptide_chain"] for r in rows}):
        subset = [r for r in rows if r["peptide_chain"] == pep]
        viable = [r for r in subset if not r["clash"]]
        exposed = [
            r for r in viable if (r["sidechain_sasa"] or 0) >= EXPOSURE_THRESHOLD
        ]
        print(f"chain {pep}: {len(viable)}/{len(subset)} rotamers clash-free")
        if viable:
            sasa = [r["sidechain_sasa"] for r in viable if r["sidechain_sasa"] is not None]
            ref_d = [r["carboxylate_to_ref"] for r in viable]
            print(f"    side-chain SASA      {min(sasa):6.1f} - {max(sasa):6.1f} A^2")
            print(f"    carboxylate -> ref   {min(ref_d):6.1f} - {max(ref_d):6.1f} A")
            print(f"    exposed (>={EXPOSURE_THRESHOLD:.0f} A^2): {len(exposed)}/{len(viable)}")
        wells = [r for r in viable if r["canonical_chi1_well"]]
        print(f"    in a canonical chi1 well: {len(wells)}")

    print(f"\nwrote {len(rows)} rows to {out.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
