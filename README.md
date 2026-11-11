# De novo minibinders against the KRAS G12D neoepitope on HLA-A\*11:01

**Anthropic × Adaptyv 2026 protein design competition — Challenge 5 (peptide–MHC). Solo entry, Track 3.**

> ### Read this first
>
> This is my **"Day 1" attempt** at creating a complete end-to-end route from
> crystal structures to 20 ranked candidate sequences. Shortcomings and potential
> improvements are noted in [`docs/shortcomings.md`](docs/shortcomings.md).
> **No design here is validated as a binder, computationally or experimentally.**
>
> **The whole project was co-executed with [Claude Code](https://claude.com/claude-code),
> with Claude Code doing the actual coding.** I set the direction, made the
> scientific decisions, chose and vetoed methods, audited structures in PyMOL,
> etc. Claude Code wrote essentially all of the Python, the tests, and the prose
> in the methods document.
>
> I intend to keep working on this — reviewing the code, revisiting the approach,
> and improving the stages one at a time. This README describes the repository as
> submitted; see [Status](#status) at the bottom.

---

## The challenge

[**Challenge page**](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/pmhc)

Design a protein that binds the KRAS G12D neoepitope `VVVGADGVGK` (residues 7–16
of mutant KRAS) presented on HLA-A\*11:01, **without** binding either the
wild-type peptide `VVVGAGGVGK` on the same HLA or the peptide-free HLA.

The objectives, hardest first:

1. **No binding to the wild-type peptide complex.**
2. **No binding to the peptide-free complex.**
3. **Strong binding to the G12D complex.**

The greatest difficulty lies in the first objective. Target and off-target differ
by one methylene-to-carboxylate substitution at a single position of a ten-residue
peptide, and the peptide contributes a minority of the complex's solvent-accessible
surface. Nearly every atom a binder could touch is identical between the two.

Affinity is the *least* important objective, and over-optimising may be actively
counterproductive — the therapeutic precedent for this exact epitope is a TCR
that discriminated cleanly at 63 µM and **acquired wild-type binding** as it was
affinity-matured to 425 nM.

## The premise

Peptide position 6 is the only chemical difference between target and off-target.
Aspartate presents a carboxylate a binder can reach; wild-type glycine presents
no side chain at all. The discriminator is therefore **aspartate versus nothing** —
presence against absence, not two chemistries to be told apart. A pocket built
to present a cation sized for that carboxylate is, against wild type, a
desolvated and unsatisfied cation: *destabilising*, not merely unrewarded. One
structural feature serves objectives 1 and 3 at once, and objective 2 is expected
to follow for free (though this ultimately depends on the epitope recognized by
the binder).

Crystallography of this system shows the mid-peptide backbone in two
arrangements: **open**, with position 6 projecting out of the groove, and
**closed**, with it turned down into the HLA F-pocket. This work designs against
the open backbone, and claims **conformational selection rather than induced
fit** — the open state need only be *populated*.

Full reasoning, including the headwinds, is in
[**`docs/methods.md`**](docs/methods.md) — the document submitted with the entry.

## What was submitted

| | |
|---|---|
| **Designs** | 20 single-chain de novo miniproteins, 61–75 aa, from 6 distinct RFdiffusion backbones |
| **Ranked by** | PRODIGY **ΔΔG = ΔG(G12D) − ΔG(wild type)**, from the same binder cofolded and minimised against peptides differing at position 6 alone |
| **Range** | −3.70 to −0.05 kcal/mol (all negative, i.e. all favour the mutant) |
| **Files** | [`submission/designs.csv`](submission/designs.csv) and [`docs/methods.md`](docs/methods.md) |

The funnel: **50 backbones → 9 passed the design gate → 90 sequences → 45 after
filtering → 87 cofolded and minimised complexes → 43 with a ΔΔG → 20 submitted.**

### The honest caveat

ΔΔG here is a differential of two PRODIGY contact-based estimates on predicted
structures. Across the 43 scored sequences the within-backbone spread (sd 1.27
kcal/mol) is *larger* than the between-backbone spread (0.69), which is what
noise looks like. More importantly, **the cofolded models bury peptide residue 6
inside the groove** (relative SASA 0.11–0.12, against 0.59 in the prepared
target), so the counter-screen could not test the interaction the designs were
built around. It was run anyway, on the argument that a buried carboxylate still
exerts electrostatic influence beyond its own solvent shell — and this system
supplies the precedent, since the natural TCR discriminates despite the residue 6
Asp pointing into the groove. Whether the ΔΔG values capture that is **not
established**. See §3.8 of the methods document.

## Pipeline

Every stage is a numbered script, driven by `make`. Run `make help` for the list.

| Stage | Command | Does | Output |
|---|---|---|---|
| 00 | `make fetch` | Download 5 PDB entries, write a provenance manifest | `data/external/` |
| 01 | `make survey` | Structural survey: open vs. closed peptide conformations | `results/survey/` |
| 02 | `make rotamers` | Gly→Asp at peptide position 6; χ1 × χ2 scan for a clash-free, exposed rotamer | `results/rotamer_scan/` |
| 03 | `make targets` | Build design targets (G12D / WT / apo, both rotamers) | `data/targets/` |
| 04 | `make relax` | Protonate at pH 7.4; restrained minimisation (amber14 + GBn2) | `data/targets_relaxed/` |
| 05 | `make nim` | **RFdiffusion** backbone generation (hosted endpoint) | `data/nim_stage05/` |
| 05 | `make backbones` | Ingest, splice onto the prepared target, apply the 5-criterion gate | `data/backbones/MANIFEST.tsv` |
| 05c | `make rotamer-sensitivity` | Re-score all 50 against the alternative Asp6 rotamer | `results/rotamer_sensitivity.tsv` |
| 06 | `make sequences` | **ProteinMPNN** sequence design, target held as fixed context | `data/sequences/` |
| 06b | `make select` | Drop low-complexity and lone-cysteine sequences | `SELECTED.tsv` |
| 07 | `make cofold` | **Boltz-2** cofolding against G12D and WT peptides | `data/cofolding_nim/` |
| 10 | `make minimise-cofold` | Protonate + restrained minimisation of each predicted complex | `*/model_min.pdb` |
| 08 | `make score-cofold-min` | PRODIGY ΔG, Levy interface classification, per-pair ipTM | `results/cofolding/SCORES_min.tsv` |
| 09 | `make package` | Rank by ΔΔG, write the submission bundle | `submission/` |

Stage numbering reflects the order stages were *planned*; 10 runs before 08 in
the final route because minimisation was added after scoring already existed.

### Definition of interface residues

Interfaces are identified by the
[Levy (2010)](https://doi.org/10.1016/j.jmb.2010.09.028) relative-SASA
classification: per-residue SASA computed in the complex and with partner atoms
removed at identical coordinates, divided by the
[Tien et al. (2013)](https://doi.org/10.1371/journal.pone.0080635) *theoretical*
maximum, classified at rSASA 0.25 into interior / surface / support / rim / core,
with **interface = support + rim + core**. Closest-approach distances are
reported descriptively and decide nothing.

It is applied over two *groups* of chains rather than chain-by-chain, because on
a three-chain system the per-chain reading counts the HLA–peptide interface and
swamps the binder signal.

The classifier is pinned to a known answer by a test on a real structure: on the
G12D receptor-bound complex it places Asp6 in **core** against the HLA
(0.572 → 0.089, 93.1 Å² lost) and in **rim** against the TCR (0.572 → 0.488,
16.2 Å²). That recovers the published indirect-electrostatic mechanism for this
system using no distances at all.

### Design gate

Five criteria, all required; **9 of 50 backbones passed.**

| # | Criterion | Threshold |
|---|---|---|
| 1 | No backbone clash | ≥ 2.0 Å to target **backbone** |
| 2 | Asp6 is an interface residue | support, rim or core |
| 3 | Asp6 buried surface area | ≥ 10 Å² |
| 4 | Peptide share of target-side interface | ≥ 10% |
| 5 | Binder compactness | Rg ≤ 1.7 × globular |

Criteria 1 and 5 both exist because of mistakes caught by inspection:
an all-atom clash rule was rejecting designs whose only close contacts were
flexible arginines, and 11 of the 20 designs that passed criteria 1–4 turned out
to be **single extended α-helices**, 88–106 Å end to end — which bury plenty of
surface while being unfoldable as standalone proteins. Both are described in
§3.6 of the methods document.

## Repository layout

```
src/kras_pmhc/      library code — structures, geometry, interface (Levy), relax,
                    build, survey, backbones, nim/mpnn clients
scripts/            one numbered script per pipeline stage
tests/              153 tests, including the 7OW6 Asp6 literature anchor
data/               targets, backbones, sequences, cofolded complexes
results/            summary tables (per-residue interface tables, scores)
submission/         the submitted CSV, FASTA, structures and methods copy
docs/
  methods.md        the submitted methods paper
  shortcomings.md   stated limitations and potential improvements
  design_plan.md    the working plan, including decisions not taken
  colab_runbook.md  the manual GPU steps, exactly as performed
notebooks/          RFdiffusion and Boltz-2 notebooks (see note below)
```

## Reproducing

```bash
conda env create -f environment.yml      # or environment.lock.yml for exact builds
conda activate kras-pmhc
make install
make test && make lint
make fetch survey rotamers targets relax
```

Stages 00–04 are fully reproducible locally. **Stages 05 and 07 are not**, for
two separate reasons:

- They need a GPU, so they ran against hosted inference endpoints. `make nim-key`
  stores an NVIDIA API key via a hidden prompt into `.nim_api_key` (mode 0600,
  gitignored — no credential is in this repository or its history).
- **The ProteinMPNN endpoint does not honour its random seed.** Two identical
  requests returned different sequences. Stage 06 therefore cannot be reproduced
  by re-running it, and the raw multi-FASTA responses are committed as the only
  record of what was designed.

`notebooks/` holds the Colab and Kaggle routes for stages 05 and 07. They are
committed for the record but the Kaggle and Boltz-2 notebooks were **never
successfully run** — free-tier GPU quota, then a Python 3.13 / DGL dependency
break documented in [`docs/colab_runbook.md`](docs/colab_runbook.md). The hosted
endpoints are what actually produced the results.

## Shortcomings

The full list, with potential improvements, is
[`docs/shortcomings.md`](docs/shortcomings.md). The four that matter most:

1. **Designs were not validated.** No refolding self-consistency check and no
   apo-condition prediction — so objective 2 was never tested at all.
2. **The aspartate is modelled, not observed.** The free complex used as the
   template carries wild-type glycine at position 6. Changing the constructed
   rotamer moves **9 of 50 gate verdicts**, with individual swings from 36.8 Å²
   of burial to 4.0 Å².
3. **All backbones were sampled at the generator's default noise**, because the
   hosted endpoint exposes no noise-scale field. The RFdiffusion authors
   recommend 0.5 for protein–protein interfaces *precisely because it improves
   designability*.
4. **Scale.** 50 backbones is small for de novo binder design, where published
   hit rates make hundreds to thousands the norm. And every surviving binder is
   84–97% helix with zero strand content.

## Status

Submitted to the competition on 9 November 2026; the repository state at
submission is tagged [`submission-day1`](../../releases/tag/submission-day1).
Work continuing beyond that tag is an improvement pass, not the entry.

## Key references

Full list in [`docs/methods.md`](docs/methods.md).

- Levy (2010) *J Mol Biol* **403**:660 — interface classification
- Tien et al. (2013) *PLoS ONE* **8**:e80635 — theoretical maximum ASA
- Watson et al. (2023) *Nature* **620**:1089 — RFdiffusion
- Dauparas et al. (2022) *Science* **378**:49 — ProteinMPNN
- Passaro et al. (2025) — Boltz-2
- Xue et al. (2016) *Bioinformatics* **32**:3676 — PRODIGY
- Poole et al. (2022) *Nat Commun* **13**:5333 — the KRAS G12D TCR precedent
