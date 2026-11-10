#!/usr/bin/env python
"""Stage 10 - protonate and minimise the cofolded complexes.

Boltz-2 returns heavy atoms at predicted positions; nothing in them has felt a
force field. PRODIGY counts atomic contacts across the interface, so its answer
depends on side-chain placement, and unminimised predicted side chains are the
weakest part of that input. This stage adds hydrogens and runs a restrained
minimisation before PRODIGY sees the structure.

**Protocol, identical to stage 04 and for the same reason.** Backbone heavy
atoms are restrained; every side chain is free. The predicted pose is the thing
under test and must not drift, while the side chains - what PRODIGY actually
counts - are free to settle. amber14 with GBn2 implicit solvent, no dynamics.

The G12D and wild-type complexes of a given sequence get a byte-identical
protocol. That is the whole basis of the ddG comparison: if the two were
relaxed differently they would differ for reasons unrelated to position 6.

Writes, per job, alongside the prediction:

    <job_id>/model_min.pdb     minimised, with hydrogens
    MINIMISATION.tsv           energies before and after, and the drift

Usage:
    python scripts/10_minimize_cofolded.py
    python scripts/10_minimize_cofolded.py --limit 2        # try a couple first
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import relax  # noqa: E402

DEFAULT_ROOT = ROOT / "data" / "cofolding_nim"
RESULTS = "RESULTS.tsv"
LOG = "MINIMISATION.tsv"

#: Above this, minimisation has moved the predicted pose rather than tidied it,
#: and the structure is reported but flagged. Backbone atoms are restrained at
#: 10 kcal/mol/A^2, so a few hundredths of an angstrom is expected.
MAX_BACKBONE_DRIFT_A = 0.5

FIELDS = (
    "job_id", "n_atoms", "n_restrained",
    "energy_before", "energy_after", "backbone_drift_a", "drifted",
    "elapsed_s",
)


def cif_to_pdb(cif_path: Path, pdb_path: Path):
    """Boltz mmCIF to PDB, which is what PDBFixer reads."""
    array = get_structure(CIFFile.read(str(cif_path)), model=1)
    out = PDBFile()
    out.set_structure(array)
    out.write(str(pdb_path))
    return array


def backbone_drift(before, after_path: Path) -> float:
    """RMSD over backbone heavy atoms, which the restraint was holding."""
    minimised = PDBFile.read(str(after_path)).get_structure(model=1)
    mask = np.isin(before.atom_name, list(relax.BACKBONE_ATOMS))
    keys = {
        (str(c), int(r), str(a)): i
        for i, (c, r, a) in enumerate(
            zip(
                before.chain_id[mask],
                before.res_id[mask],
                before.atom_name[mask],
                strict=True,
            )
        )
    }
    before_coords = before.coord[mask]

    mask2 = np.isin(minimised.atom_name, list(relax.BACKBONE_ATOMS))
    after_coords, matched = [], []
    for coord, chain, res, atom in zip(
        minimised.coord[mask2],
        minimised.chain_id[mask2],
        minimised.res_id[mask2],
        minimised.atom_name[mask2],
        strict=True,
    ):
        key = (str(chain), int(res), str(atom))
        if key in keys:
            after_coords.append(coord)
            matched.append(keys[key])
    if not matched:
        return float("nan")
    a = np.asarray(after_coords)
    b = before_coords[matched]
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def read_log(root: Path) -> set[str]:
    path = root / LOG
    if not path.is_file():
        return set()
    with path.open() as handle:
        return {row["job_id"] for row in csv.DictReader(handle, delimiter="\t")}


def append_log(root: Path, row: dict) -> None:
    path = root / LOG
    exists = path.is_file()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--limit", type=int, help="minimise at most this many")
    parser.add_argument("--platform", default="CPU")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    results = args.root / RESULTS
    if not results.is_file():
        raise SystemExit(f"{results} not found")
    with results.open() as handle:
        jobs = list(csv.DictReader(handle, delimiter="\t"))

    done = set() if args.force else read_log(args.root)
    pending = [
        job for job in jobs
        if (args.root / job["job_id"] / "model.cif").is_file()
        and (args.force or job["job_id"] not in done)
    ]
    if args.limit:
        pending = pending[: args.limit]

    print(f"{len(jobs)} prediction(s), {len(done)} already minimised, "
          f"{len(pending)} to do")
    print(f"amber14 + GBn2, backbone restrained at "
          f"{relax.DEFAULT_RESTRAINT_K} kcal/mol/A^2, platform {args.platform}\n")

    drifted = []
    for index, job in enumerate(pending, start=1):
        job_dir = args.root / job["job_id"]
        started = time.time()

        original = cif_to_pdb(job_dir / "model.cif", job_dir / "model_heavy.pdb")
        topology, protonated = relax.protonate(job_dir / "model_heavy.pdb")
        positions, before, after, n_restrained = relax.minimise(
            topology, protonated, platform_name=args.platform
        )
        out = job_dir / "model_min.pdb"
        relax.write_pdb(topology, positions, out)

        drift = backbone_drift(original, out)
        has_drifted = bool(drift > MAX_BACKBONE_DRIFT_A)
        if has_drifted:
            drifted.append(job["job_id"])
        elapsed = time.time() - started

        append_log(args.root, {
            "job_id": job["job_id"],
            "n_atoms": topology.getNumAtoms(),
            "n_restrained": n_restrained,
            "energy_before": round(before, 1),
            "energy_after": round(after, 1),
            "backbone_drift_a": round(drift, 4),
            "drifted": has_drifted,
            "elapsed_s": round(elapsed, 1),
        })
        print(f"[{index}/{len(pending)}] {job['job_id']:<22} "
              f"{before:>12.0f} -> {after:>10.0f} kcal/mol  "
              f"drift {drift:.3f} A  {elapsed:.0f}s"
              + ("  DRIFTED" if has_drifted else ""), flush=True)

    if drifted:
        print(f"\n{len(drifted)} structure(s) moved more than "
              f"{MAX_BACKBONE_DRIFT_A} A despite the restraint: {drifted}")
    print(f"\nminimised structures: {args.root}/<job_id>/model_min.pdb")
    print(f"log: {args.root / LOG}")
    print("\nScore them with:")
    print("  python scripts/08_score_cofolded.py --structure model_min.pdb "
          "--out results/cofolding/SCORES_min.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
