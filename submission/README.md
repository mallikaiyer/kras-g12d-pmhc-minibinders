# KRAS G12D / HLA-A*11:01 minibinders

10 designed binder sequences from 7 distinct backbones.
Generated 2026-11-09 23:02 UTC. Full method in `docs/methods.md`; pipeline in this repository.

## Target

| | |
|---|---|
| Peptide | `VVVGADGVGK` (KRAS G12D, residues 7-16) |
| MHC class I | HLA-A*11:01 heavy chain, 281 aa |
| Light chain | beta-2 microglobulin, 99 aa |
| Designed against | heavy-chain residues 1-180 (alpha1/alpha2) + peptide |

The wild-type peptide is `VVVGAGGVGK`; it differs from the target at position
6 alone, which is the residue every design is aimed at.

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
