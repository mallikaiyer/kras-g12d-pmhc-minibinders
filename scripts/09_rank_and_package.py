#!/usr/bin/env python
"""Stage 09 - rank the designs and write the submission package.

The ranking is lexicographic, in the order of precedence set for this project:

1. **Mutation specificity** - how much of Asp6 the binder buries.
2. **Peptide-dependent binding** - the peptide's share of the target-side
   interface.
3. **Affinity** - last, and only among designs that clear the first two.

All three are measured on the **crystal-derived** design target from stage 04,
not on the cofolded models. That is a deliberate choice and the reason is
recorded in docs/methods.md: the cofolded complexes bury Asp6 inside the groove
(rSASA 0.11-0.12 against 0.59 in the crystal-derived target), while 7OW6 shows a
TCR contacting Asp6 directly. Ranking on structures that do not present the
residue under test would rank on an artefact. The cofolding numbers are carried
into the table as evidence, flagged, and not ranked on.

Sequence choice within a design is by ProteinMPNN score, and at most two
sequences per design are submitted so that the set is not dominated by one
backbone.

Outputs:
    submission/designs.fasta   ranked binder sequences
    submission/designs.tsv     the evidence behind every row
    submission/README.md       what the package is and how it was made

Usage:
    python scripts/09_rank_and_package.py
    python scripts/09_rank_and_package.py --max-designs 10 --per-design 2
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

from kras_pmhc.constants import (  # noqa: E402
    B2M_SEQUENCE,
    HLA_A1101_HEAVY_CHAIN,
    MUTANT_POSITION,
)

MANIFEST = ROOT / "data" / "backbones" / "MANIFEST.tsv"
SELECTED = ROOT / "data" / "sequences" / "SELECTED.tsv"
ROTAMER = ROOT / "results" / "rotamer_sensitivity.tsv"
COFOLD = ROOT / "results" / "cofolding" / "SCORES.tsv"
OUT_DIR = ROOT / "submission"

PEPTIDE_G12D = "VVVGADGVGK"

FIELDS = (
    "rank", "design_id", "backbone_id", "sample", "binder_length",
    "binder_sequence",
    "asp6_class", "asp6_sasa_loss", "peptide_share", "compactness",
    "mpnn_score", "rotamer_robust",
    "cofold_iptm_binder_peptide", "cofold_prodigy_dg", "cofold_caveat",
)


def read_tsv(path: Path, required: bool = True) -> list[dict]:
    if not path.is_file():
        if required:
            raise SystemExit(f"{path} not found")
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-designs", type=int, default=10)
    parser.add_argument("--per-design", type=int, default=2)
    args = parser.parse_args()

    backbones = {
        row["backbone_id"]: row
        for row in read_tsv(MANIFEST)
        if row["passes_gate"] == "True"
    }
    sequences = [
        row for row in read_tsv(SELECTED) if row["backbone_id"] in backbones
    ]
    if not sequences:
        raise SystemExit("no sequence survives the gate and the selection filters")

    rotamer = {row["backbone_id"]: row for row in read_tsv(ROTAMER, required=False)}
    cofold = collections.defaultdict(dict)
    for row in read_tsv(COFOLD, required=False):
        cofold[(row["backbone_id"], row["sample"])][row["condition"]] = row

    # 1 and 2: order the backbones. Asp6 burial first, peptide share second.
    ranked_backbones = sorted(
        backbones.values(),
        key=lambda r: (-float(r["asp6_sasa_loss"]), -float(r["peptide_share"])),
    )

    by_design = collections.defaultdict(list)
    for row in sequences:
        by_design[row["backbone_id"]].append(row)

    # 3: affinity last, and only as the tie-break within a backbone.
    candidates = {
        backbone["backbone_id"]: sorted(
            by_design.get(backbone["backbone_id"], []),
            key=lambda r: float(r["score"]),
        )[: args.per_design]
        for backbone in ranked_backbones
    }

    # Allocated round-robin rather than by filling each backbone in turn: with
    # 10 slots and 2 per design, filling in order would spend every slot on the
    # top five backbones and leave two unrepresented. One sequence from each
    # backbone comes first, then second sequences in the same order. Submitting
    # several sequences from one backbone is a bet on that backbone, and there
    # is no evidence here strong enough to justify concentrating the set.
    order = []
    for slot in range(args.per_design):
        for backbone in ranked_backbones:
            picks = candidates.get(backbone["backbone_id"], [])
            if slot < len(picks):
                order.append((backbone, picks[slot]))

    rows, rank = [], 0
    for backbone, candidate in order:
        design_id = backbone["backbone_id"]
        if rank >= args.max_designs:
            break
        rank += 1
        g12d = cofold.get((design_id, candidate["sample"]), {}).get("g12d", {})
        rows.append({
            "rank": rank,
            "design_id": f"kras_g12d_binder_{rank:02d}",
            "backbone_id": design_id,
            "sample": candidate["sample"],
            "binder_length": len(candidate["sequence"]),
            "binder_sequence": candidate["sequence"],
            "asp6_class": backbone["asp6_class"],
            "asp6_sasa_loss": backbone["asp6_sasa_loss"],
            "peptide_share": backbone["peptide_share"],
            "compactness": backbone["compactness"],
            "mpnn_score": candidate["score"],
            "rotamer_robust": rotamer.get(design_id, {}).get(
                "alternative_passes", "unknown"
            ),
            "cofold_iptm_binder_peptide": g12d.get("iptm_binder_peptide", ""),
            "cofold_prodigy_dg": g12d.get("prodigy_dg", ""),
            "cofold_caveat": (
                "Asp6 buried in cofolded model; not ranked on" if g12d else "not cofolded"
            ),
        })

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with (OUT_DIR / "designs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    fasta = []
    for row in rows:
        fasta.append(
            f">{row['design_id']} backbone={row['backbone_id']} "
            f"mpnn_sample={row['sample']} length={row['binder_length']} "
            f"asp6_class={row['asp6_class']} "
            f"asp6_sasa_loss={row['asp6_sasa_loss']}A2 "
            f"peptide_share={row['peptide_share']}"
        )
        fasta.append(row["binder_sequence"])
    (OUT_DIR / "designs.fasta").write_text("\n".join(fasta) + "\n")

    stamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    designs = len({row["backbone_id"] for row in rows})
    readme = f"""# KRAS G12D / HLA-A*11:01 minibinders

{len(rows)} designed binder sequences from {designs} distinct backbones.
Generated {stamp}. Full method in `docs/methods.md`; pipeline in this repository.

## Target

| | |
|---|---|
| Peptide | `{PEPTIDE_G12D}` (KRAS G12D, residues 7-16) |
| MHC class I | HLA-A*11:01 heavy chain, {len(HLA_A1101_HEAVY_CHAIN)} aa |
| Light chain | beta-2 microglobulin, {len(B2M_SEQUENCE)} aa |
| Designed against | heavy-chain residues 1-180 (alpha1/alpha2) + peptide |

The wild-type peptide is `VVVGAGGVGK`; it differs from the target at position
{MUTANT_POSITION} alone, which is the residue every design is aimed at.

## Files

| file | contents |
|---|---|
| `designs.fasta` | the binder sequences, ranked |
| `designs.tsv` | the evidence behind each row |

## How the ranking was made

Lexicographic, in this order:

1. **Mutation specificity** - Asp6 buried surface area on binding, by the Levy
   (2010) rSASA classification. All designs bury >= 10 A2.
2. **Peptide-dependent binding** - the peptide's share of target-side interface
   area, >= 10% for every design here.
3. **Affinity** - ProteinMPNN score, used only to choose between sequences on
   the same backbone.

## What these numbers are not

Read `docs/methods.md` before treating any of this as validated. In short:

* Designs are **not** experimentally or computationally validated as binders.
  Cofolding placed the binder interface at ipTM 0.16-0.54 while the same
  predictions put the real HLA-peptide interface at 0.94-0.96.
* The cofolded models **bury Asp6** inside the groove (rSASA 0.11-0.12) where
  the crystal-derived target has it exposed (0.59), so that screen could not
  test the intended interaction and is reported but not ranked on.
* Interface numbers are measured against a **constructed** Asp6 side chain;
  7OW3 carries wild-type Gly at that position. Changing the rotamer moves 9 of
  50 gate verdicts, and `rotamer_robust` records which designs survive the
  alternative.
* All backbones were sampled at RFdiffusion's **default noise** because the
  hosted endpoint used exposes no noise-scale field; the authors recommend 0.5
  for interfaces precisely because it improves designability.
"""
    (OUT_DIR / "README.md").write_text(readme)

    print(f"{len(rows)} sequences from {designs} backbones\n")
    print(f"{'rank':<5}{'design':<28}{'backbone':<10}{'Asp6 loss':>10}{'share':>7}{'mpnn':>7}")
    print("=" * 67)
    for row in rows:
        print(f"{row['rank']:<5}{row['design_id']:<28}{row['backbone_id']:<10}"
              f"{float(row['asp6_sasa_loss']):>10.1f}"
              f"{float(row['peptide_share']):>7.0%}{float(row['mpnn_score']):>7.3f}")
    print(f"\nwrote {(OUT_DIR / 'designs.fasta').relative_to(ROOT)}")
    print(f"      {(OUT_DIR / 'designs.tsv').relative_to(ROOT)}")
    print(f"      {(OUT_DIR / 'README.md').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
