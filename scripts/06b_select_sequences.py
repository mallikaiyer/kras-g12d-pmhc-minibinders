#!/usr/bin/env python
"""Stage 06b - select the designed sequences worth carrying forward.

Two exclusions, both decided by the user on 9 Nov 2026:

**Low complexity.** A sequence with one residue type at 40% or more. These are
not a sampling artefact to be averaged away - they are what ProteinMPNN returns
when a backbone has no core to pack, and they were the first evidence that many
stage-05 designs were single extended helices. nim3_0's best-scoring sequence is
87% one residue and 88% hydrophobic.

**A sole cysteine.** Every cysteine-containing sequence here has exactly one, so
none can form an intramolecular disulfide. An unpaired cysteine in a binder
destined for expression can scramble with host thiols or drive aggregation, and
it buys nothing in a design that was never meant to have a disulfide.

SEQUENCES.tsv is left intact as the full record of what was designed; this
writes the surviving subset alongside it. The split is kept rather than filtered
in place because the discarded sequences are evidence about the backbones.

Usage:
    python scripts/06b_select_sequences.py
"""

from __future__ import annotations

import argparse
import collections
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEQUENCES = ROOT / "data" / "sequences" / "SEQUENCES.tsv"
SELECTED = ROOT / "data" / "sequences" / "SELECTED.tsv"

#: One residue type at or above this fraction of the sequence.
MAX_SINGLE_RESIDUE_FRACTION = 0.40

#: Cysteines tolerated. Zero: a lone cysteine is a liability and these
#: sequences have no paired ones.
MAX_CYSTEINES = 0


def load() -> list[dict]:
    if not SEQUENCES.is_file():
        raise SystemExit(f"{SEQUENCES} not found - run scripts/06_mpnn_sequences.py")
    with SEQUENCES.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def reasons(row: dict) -> list[str]:
    """Why this sequence is excluded, if it is."""
    out = []
    if float(row["most_common_fraction"]) >= MAX_SINGLE_RESIDUE_FRACTION:
        out.append(
            f"low complexity ({row['most_common_residue']} "
            f"{float(row['most_common_fraction']):.0%})"
        )
    if int(row["n_cysteine"]) > MAX_CYSTEINES:
        out.append(f"{row['n_cysteine']} cysteine")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.parse_args()

    rows = load()
    kept, dropped = [], []
    for row in rows:
        why = reasons(row)
        (dropped if why else kept).append((row, why))

    if not kept:
        raise SystemExit("every sequence was excluded; nothing to carry forward")

    with SELECTED.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(row for row, _ in kept)

    by_design = collections.Counter(row["backbone_id"] for row, _ in kept)
    all_designs = {row["backbone_id"] for row in rows}
    lost = sorted(all_designs - set(by_design))

    print(f"kept {len(kept)} of {len(rows)} sequences, across "
          f"{len(by_design)} of {len(all_designs)} designs\n")
    for design, n in sorted(by_design.items()):
        print(f"  {design:<9} {n:>2}")
    if lost:
        print(f"\nno usable sequence at all: {', '.join(lost)}")

    counts = collections.Counter()
    for _, why in dropped:
        for reason in why:
            counts["low complexity" if reason.startswith("low") else "cysteine"] += 1
    print(f"\nexcluded {len(dropped)} sequence(s): "
          + ", ".join(f"{n} for {reason}" for reason, n in counts.most_common()))
    print(f"\nwrote {SELECTED.relative_to(ROOT)}")
    print(f"{SEQUENCES.name} is unchanged: the excluded sequences are evidence "
          "about their backbones, not noise")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
