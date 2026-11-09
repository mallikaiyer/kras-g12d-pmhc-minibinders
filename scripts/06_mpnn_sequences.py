#!/usr/bin/env python
"""Stage 06 - design sequences onto the backbones that passed the stage-05 gate.

One API call per backbone, ten sequences each. The binder chain is designed and
the HLA and peptide are held as fixed context, so ProteinMPNN is conditioned on
the target's sequence as well as its shape.

The raw multi-FASTA response is written verbatim alongside the parsed table,
because the hosted endpoint ignores ``random_seed``: re-running this script
produces different sequences, so the response is the record. Nothing here
overwrites an existing response unless --force is given.

Outputs:
    results/sequences/<backbone_id>.fasta   raw response, exactly as returned
    data/sequences/SEQUENCES.tsv            one row per designed sequence

Usage:
    python scripts/06_mpnn_sequences.py
    python scripts/06_mpnn_sequences.py --num-sequences 10 --temperature 0.1
"""

from __future__ import annotations

import argparse
import collections
import csv
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import mpnn, nim  # noqa: E402

MANIFEST = ROOT / "data" / "backbones" / "MANIFEST.tsv"
BACKBONE_DIR = ROOT / "data" / "backbones"
FASTA_DIR = ROOT / "results" / "sequences"
OUT = ROOT / "data" / "sequences" / "SEQUENCES.tsv"
KEY_FILE = ROOT / ".nim_api_key"

#: The generated chain. The target keeps A and C and is held fixed.
BINDER_CHAIN = "B"

#: Residues flagged in the output for the user's attention rather than filtered.
#: A free cysteine is a liability in a binder destined for expression - it can
#: scramble disulfides or drive aggregation - but whether to exclude it is a
#: design decision, so it is counted and reported, not silently omitted.
WATCH_RESIDUES = ("C",)

#: Low-complexity threshold for reporting. A single residue type over this
#: fraction is the poly-alanine signature that exposed the extended-helix
#: problem at stage 05, and is worth seeing per sequence.
LOW_COMPLEXITY_FRACTION = 0.40

HYDROPHOBIC = set("AVLIMFWYC")

FIELDS = (
    "backbone_id", "sample", "temperature", "length", "sequence",
    "score", "global_score", "seq_recovery",
    "hydrophobic_fraction", "most_common_residue", "most_common_fraction",
    "n_cysteine", "low_complexity", "designed_utc",
)


def read_manifest() -> list[dict]:
    if not MANIFEST.is_file():
        raise SystemExit(f"{MANIFEST} not found - run scripts/05_ingest_backbones.py")
    with MANIFEST.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    passing = [row for row in rows if row["passes_gate"] == "True"]
    if not passing:
        raise SystemExit("no backbone passed the stage-05 gate; nothing to design")
    return passing


def profile(sequence: str) -> dict:
    counts = collections.Counter(sequence)
    residue, count = counts.most_common(1)[0]
    fraction = count / len(sequence)
    return {
        "hydrophobic_fraction": round(
            sum(counts[a] for a in HYDROPHOBIC) / len(sequence), 3
        ),
        "most_common_residue": residue,
        "most_common_fraction": round(fraction, 3),
        "n_cysteine": sum(counts[a] for a in WATCH_RESIDUES),
        "low_complexity": fraction >= LOW_COMPLEXITY_FRACTION,
    }


def append(rows: list[dict]) -> None:
    """Append after every backbone, so an interruption keeps what was designed."""
    exists = OUT.is_file()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def already_designed() -> set[str]:
    if not OUT.is_file():
        return set()
    with OUT.open() as handle:
        return {row["backbone_id"] for row in csv.DictReader(handle, delimiter="\t")}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--num-sequences", type=int, default=mpnn.DEFAULT_NUM_SEQUENCES)
    parser.add_argument("--temperature", type=float, default=mpnn.DEFAULT_TEMPERATURE)
    parser.add_argument(
        "--force",
        action="store_true",
        help="redesign backbones already present in SEQUENCES.tsv. The service "
             "ignores the seed, so this yields different sequences and discards "
             "the only record of the previous ones",
    )
    args = parser.parse_args()

    key = nim.read_api_key(KEY_FILE)
    passing = read_manifest()
    done = set() if args.force else already_designed()
    FASTA_DIR.mkdir(parents=True, exist_ok=True)

    print(f"{len(passing)} backbone(s) passed the gate")
    print(f"{args.num_sequences} sequences each at T={args.temperature}, "
          f"designing chain {BINDER_CHAIN} with the target fixed")
    if done:
        print(f"skipping {len(done)} already designed (use --force to redo)")

    total = 0
    for row in passing:
        backbone_id = row["backbone_id"]
        if backbone_id in done:
            continue

        pdb = BACKBONE_DIR / f"{backbone_id}.pdb"
        expected = int(row["binder_length"])
        print(f"\n{backbone_id}  ({expected} aa)")

        try:
            designs, mfasta = mpnn.design(
                key,
                pdb.read_text(),
                [BINDER_CHAIN],
                num_sequences=args.num_sequences,
                temperature=args.temperature,
                expected_length=expected,
            )
        except (mpnn.MpnnError, nim.NimError) as error:
            print(f"  FAILED: {error}")
            continue

        (FASTA_DIR / f"{backbone_id}.fasta").write_text(mfasta)

        stamp = datetime.now(UTC).isoformat(timespec="seconds")
        rows = []
        for design in designs:
            rows.append({
                "backbone_id": backbone_id,
                "sample": design.sample,
                "temperature": design.temperature,
                "length": design.length,
                "sequence": design.sequence,
                "score": round(design.score, 4),
                "global_score": round(design.global_score, 4),
                "seq_recovery": round(design.seq_recovery, 4),
                "designed_utc": stamp,
                **profile(design.sequence),
            })
        append(rows)
        total += len(rows)

        best = min(designs, key=lambda d: d.score)
        flags = [r for r in rows if r["low_complexity"]]
        cys = sum(r["n_cysteine"] for r in rows)
        print(f"  {len(designs)} sequences, best score {best.score:.3f}")
        print(f"  {best.sequence}")
        if flags:
            print(f"  {len(flags)} low-complexity (>= "
                  f"{LOW_COMPLEXITY_FRACTION:.0%} one residue)")
        if cys:
            print(f"  {cys} cysteine(s) across the set - decide before submission")

    print(f"\n{total} sequences written to {OUT.relative_to(ROOT)}")
    print(f"raw responses in {FASTA_DIR.relative_to(ROOT)}/")
    return 0 if total else 1


if __name__ == "__main__":
    raise SystemExit(main())
