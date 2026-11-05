#!/usr/bin/env python
"""Stage 00 - fetch the five relevant PDB entries from RCSB.

Downloads are kept out of git (they are regenerable from the accession codes)
but a manifest with sizes and SHA-256 digests is committed, so that a later
re-run can be shown to have analysed identical files.

Usage:
    python scripts/00_fetch_structures.py [--force]
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

from biotite.database import rcsb

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from kras_pmhc.constants import STRUCTURES  # noqa: E402

PDB_DIR = REPO_ROOT / "data" / "external" / "pdb"
MANIFEST = REPO_ROOT / "data" / "external" / "MANIFEST.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--force", action="store_true", help="re-download entries already present"
    )
    args = parser.parse_args()

    PDB_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    for pdb_id, meta in sorted(STRUCTURES.items()):
        target = PDB_DIR / f"{pdb_id.lower()}.cif"
        if target.exists() and not args.force:
            print(f"  have  {pdb_id}  {meta.description}")
        else:
            print(f"  fetch {pdb_id}  {meta.description}")
            # biotite writes <pdb_id>.cif into the directory, lower-cased.
            rcsb.fetch(pdb_id, "cif", target_path=str(PDB_DIR), overwrite=True)
            fetched = PDB_DIR / f"{pdb_id}.cif"
            if fetched.exists() and fetched != target:
                fetched.rename(target)

        if not target.exists():
            print(f"  ERROR {pdb_id}: expected {target} after fetch", file=sys.stderr)
            return 1

        rows.append(
            {
                "pdb_id": pdb_id,
                "filename": target.name,
                "bytes": target.stat().st_size,
                "sha256": sha256(target),
                "peptide": meta.peptide,
                "has_tcr": meta.has_tcr,
                "description": meta.description,
            }
        )

    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n{len(rows)} entries in {PDB_DIR}")
    print(f"manifest written to {MANIFEST.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
