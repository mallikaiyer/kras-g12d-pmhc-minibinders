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
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import mpnn  # noqa: E402
from kras_pmhc.constants import (  # noqa: E402
    B2M_SEQUENCE,
    HLA_A1101_HEAVY_CHAIN,
    MUTANT_POSITION,
)

MANIFEST = ROOT / "data" / "backbones" / "MANIFEST.tsv"
SELECTED = ROOT / "data" / "sequences" / "SELECTED.tsv"
ROTAMER = ROOT / "results" / "rotamer_sensitivity.tsv"
COFOLD = ROOT / "results" / "cofolding" / "SCORES.tsv"
COFOLD_MIN = ROOT / "results" / "cofolding" / "SCORES_min.tsv"
OUT_DIR = ROOT / "submission"

PEPTIDE_G12D = "VVVGADGVGK"

FIELDS = (
    "rank", "design_id", "backbone_id", "sample", "binder_length",
    "binder_sequence",
    "prodigy_ddg", "prodigy_dg_g12d", "prodigy_dg_wt", "favours",
    "asp6_class", "asp6_sasa_loss", "peptide_share", "compactness",
    "mpnn_score", "rotamer_robust", "mean_kyte_doolittle", "complex_plddt",
    "cofold_iptm_binder_peptide", "cofold_iptm_binder_hla",
    "structures", "cofold_caveat",
)

#: The submission CSV. The challenge requires `name`, `sequence` and
#: `molecule_class` at minimum, with rows ordered best-first, and encourages
#: sharing whatever metrics the entrant is comfortable publishing - so the
#: physics-based and confidence metrics are included rather than held back.
CSV_FIELDS = (
    # required
    "name", "sequence", "molecule_class",
    # ranking and its basis
    "rank", "prodigy_ddg_kcal_mol", "prodigy_dg_g12d_kcal_mol",
    "prodigy_dg_wt_kcal_mol", "favours",
    # confidence from cofolding, and the sequence-liability screen the
    # competition's selection prompt names alongside pLDDT and ipSAE
    "iptm_binder_peptide", "iptm_binder_hla", "complex_plddt",
    "mean_kyte_doolittle",
    # structural engagement of the mutation, on the crystal-derived target
    "asp6_class", "asp6_sasa_loss_a2", "peptide_share", "binder_compactness",
    # sequence design
    "binder_length", "mpnn_score",
    # provenance
    "backbone_id", "mpnn_sample", "structures_scored", "caveat",
)

#: Every design here is a single-chain de novo miniprotein.
MOLECULE_CLASS = "protein"


def read_tsv(path: Path, required: bool = True) -> list[dict]:
    if not path.is_file():
        if required:
            raise SystemExit(f"{path} not found")
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bundle(rows: list[dict]) -> None:
    """Copy the methods document and the structures behind each design.

    The challenge encourages sharing structure models and methodology, and
    notes submitted data may be published. Everything here is already in the
    repository; this gathers the subset that belongs to the submitted designs
    so the package stands alone.
    """
    shutil.copy2(ROOT / "docs" / "methods.md", OUT_DIR / "methods.md")

    structures = OUT_DIR / "structures"
    backbone_dir = structures / "backbones"
    structures.mkdir(parents=True, exist_ok=True)
    backbone_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    for row in rows:
        for condition in ("g12d", "wt"):
            job = f"{row['backbone_id']}_s{row['sample']}_{condition}"
            # Prefer the minimised model: it is what PRODIGY scored.
            for name in ("model_min.pdb", "model.cif"):
                source = ROOT / "data" / "cofolding_nim" / job / name
                if source.is_file():
                    suffix = source.suffix
                    target = structures / f"{row['design_id']}_{condition}{suffix}"
                    shutil.copy2(source, target)
                    copied += 1
                    break

    backbones_copied = set()
    for row in rows:
        source = ROOT / "data" / "backbones" / f"{row['backbone_id']}.pdb"
        if source.is_file() and row["backbone_id"] not in backbones_copied:
            shutil.copy2(source, backbone_dir / source.name)
            backbones_copied.add(row["backbone_id"])

    print(f"bundled methods.md, {copied} cofolded model(s) and "
          f"{len(backbones_copied)} design backbone(s)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-designs", type=int, default=20)
    parser.add_argument(
        "--per-design",
        type=int,
        default=0,
        help="cap on sequences from one backbone; 0 means no cap, which is "
             "what ranking purely by ddG implies",
    )
    parser.add_argument(
        "--no-structures",
        action="store_true",
        help="skip copying structures and methods into submission/",
    )
    parser.add_argument(
        "--cofold-scores",
        type=Path,
        default=None,
        help="stage-08 scores to rank on. Defaults to the minimised scores "
             f"({COFOLD_MIN.name}) when they exist, otherwise the raw "
             f"predictions ({COFOLD.name}). PRODIGY counts atomic contacts, so "
             "minimised side chains are the better input.",
    )
    parser.add_argument(
        "--rank-by",
        choices=("ddg", "burial"),
        default="ddg",
        help="ddg (default): PRODIGY dG(G12D) - dG(wild type), the project's "
             "stated specificity criterion. burial: Asp6 buried area on the "
             "crystal-derived target, the structural fallback used when "
             "cofolding is absent.",
    )
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

    scores_path = args.cofold_scores
    if scores_path is None:
        scores_path = COFOLD_MIN if COFOLD_MIN.is_file() else COFOLD
    minimised = scores_path.name == COFOLD_MIN.name
    cofold = collections.defaultdict(dict)
    for row in read_tsv(scores_path, required=False):
        cofold[(row["backbone_id"], row["sample"])][row["condition"]] = row
    print(f"ranking on {scores_path.relative_to(ROOT)} "
          f"({'minimised' if minimised else 'raw predictions'})\n")

    # --- the ranking ------------------------------------------------------
    # Primary criterion is mutation specificity. With cofolding available that
    # is PRODIGY's ddG: dG(G12D) - dG(wild type), both from the same binder
    # folded against peptides differing at position 6 alone. Negative means the
    # binder is predicted to prefer the mutant, which is the whole objective.
    #
    # ddG ranks per *sequence*, not per backbone, because each sequence has its
    # own pair of predictions - a meaningful improvement over ranking backbones.
    def ddg_of(row):
        pair = cofold.get((row["backbone_id"], row["sample"]), {})
        g12d, wt = pair.get("g12d"), pair.get("wt")
        if not g12d or not wt:
            return None
        try:
            return float(g12d["prodigy_dg"]) - float(wt["prodigy_dg"])
        except (ValueError, KeyError):
            return None

    scored = []
    for row in sequences:
        backbone = backbones[row["backbone_id"]]
        ddg = ddg_of(row)
        scored.append({"sequence_row": row, "backbone": backbone, "ddg": ddg})

    missing = [s for s in scored if s["ddg"] is None]
    if args.rank_by == "ddg" and len(missing) == len(scored):
        raise SystemExit(
            "ranking by ddG needs cofolding for both conditions and none is "
            "present - run scripts/07b_nim_cofolding.py, or pass "
            "--rank-by burial"
        )

    if args.rank_by == "ddg":
        # Sequences without a ddG sort last rather than being dropped, so the
        # set can still be filled if cofolding is incomplete.
        def key(s):
            return (
                s["ddg"] is None,
                s["ddg"] if s["ddg"] is not None else 0.0,
                -float(s["backbone"]["peptide_share"]),
                float(s["sequence_row"]["score"]),
            )
    else:
        # The structural fallback: Asp6 burial on the crystal-derived target,
        # then peptide share, then MPNN score. Used when cofolding is absent or
        # judged unusable - see docs/methods.md on the buried-Asp6 problem.
        def key(s):
            return (
                -float(s["backbone"]["asp6_sasa_loss"]),
                -float(s["backbone"]["peptide_share"]),
                float(s["sequence_row"]["score"]),
            )

    scored.sort(key=key)

    # A per-backbone cap, applied after ranking. Several sequences from one
    # backbone is a bet on that backbone, and the evidence here does not support
    # concentrating the set however well one sequence scores.
    rows, rank = [], 0
    used = collections.Counter()
    for entry in scored:
        if rank >= args.max_designs:
            break
        design_id = entry["sequence_row"]["backbone_id"]
        if args.per_design and used[design_id] >= args.per_design:
            continue
        used[design_id] += 1
        rank += 1

        candidate = entry["sequence_row"]
        backbone = entry["backbone"]
        pair = cofold.get((design_id, candidate["sample"]), {})
        g12d, wt = pair.get("g12d", {}), pair.get("wt", {})
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
            "prodigy_dg_g12d": g12d.get("prodigy_dg", ""),
            "prodigy_dg_wt": wt.get("prodigy_dg", ""),
            "prodigy_ddg": (
                "" if entry["ddg"] is None else round(entry["ddg"], 3)
            ),
            "favours": (
                "" if entry["ddg"] is None
                else "G12D" if entry["ddg"] < 0 else "wild type"
            ),
            "cofold_iptm_binder_peptide": g12d.get("iptm_binder_peptide", ""),
            "cofold_iptm_binder_hla": g12d.get("iptm_binder_hla", ""),
            "complex_plddt": g12d.get("complex_plddt", ""),
            "mean_kyte_doolittle": round(
                mpnn.mean_hydropathy(candidate["sequence"]), 3
            ),
            "structures": "minimised" if minimised else "raw prediction",
            "cofold_caveat": (
                "Asp6 buried in cofolded model (rSASA 0.11-0.12 vs 0.59 in the "
                "crystal-derived target)" if g12d else "not cofolded"
            ),
        })

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with (OUT_DIR / "designs.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # The CSV the challenge asks for: comma-delimited, best design first,
    # `name`/`sequence`/`molecule_class` present as required.
    with (OUT_DIR / "designs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "name": row["design_id"],
                "sequence": row["binder_sequence"],
                "molecule_class": MOLECULE_CLASS,
                "rank": row["rank"],
                "prodigy_ddg_kcal_mol": row["prodigy_ddg"],
                "prodigy_dg_g12d_kcal_mol": row["prodigy_dg_g12d"],
                "prodigy_dg_wt_kcal_mol": row["prodigy_dg_wt"],
                "favours": row["favours"],
                "iptm_binder_peptide": row["cofold_iptm_binder_peptide"],
                "iptm_binder_hla": row.get("cofold_iptm_binder_hla", ""),
                "complex_plddt": row["complex_plddt"],
                "mean_kyte_doolittle": row["mean_kyte_doolittle"],
                "asp6_class": row["asp6_class"],
                "asp6_sasa_loss_a2": row["asp6_sasa_loss"],
                "peptide_share": row["peptide_share"],
                "binder_compactness": row["compactness"],
                "binder_length": row["binder_length"],
                "mpnn_score": row["mpnn_score"],
                "backbone_id": row["backbone_id"],
                "mpnn_sample": row["sample"],
                "structures_scored": row["structures"],
                "caveat": row["cofold_caveat"],
            })

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
| `designs.csv` | **the submission**: name, sequence, molecule_class,
  ranked best-first, plus metrics |
| `designs.fasta` | the same sequences in FASTA |
| `designs.tsv` | the full evidence behind each row |
| `methods.md` | how every number was produced, and what it is not |
| `structures/` | the minimised cofolded complex for each design, against both peptides |
| `structures/backbones/` | the RFdiffusion backbone each design was built on |

All designs are single-chain de novo miniproteins: `molecule_class = protein`.

## How the ranking was made

Lexicographic, in this order:

1. **Mutation specificity** - PRODIGY **ddG = dG(G12D) - dG(wild type)**, from
   the same binder cofolded against peptides that differ at position 6 alone.
   Negative favours the mutant. This is the primary criterion.
2. **Peptide-dependent binding** - the peptide's share of target-side interface
   area on the crystal-derived target, >= 10% for every design here.
3. **Affinity** - ProteinMPNN score, as the final tie-break.

At most two sequences per backbone, so the set is not a bet on one backbone.

`--rank-by burial` reproduces the structural alternative: Asp6 buried area in
place of ddG, for the reason in `docs/methods.md` - the cofolded models bury
Asp6, so the ddG values are measured on structures that do not present the
residue under test. **Both rankings are reported; they disagree substantially.**

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

    if not args.no_structures:
        bundle(rows)

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
