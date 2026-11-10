# Methods — de novo minibinders against KRAS G12D on HLA-A\*11:01

Anthropic × Adaptyv 2026, Challenge 5. Solo entry, Track 3.

Everything below is reproducible from this repository: each stage is a numbered
script under `scripts/`, driven by `make`, with 149 tests covering the analysis
code. Where a step could not be scripted — the GPU notebooks — the manual
actions are recorded in `docs/colab_runbook.md`.

---

## 1. Target and premise

| | |
|---|---|
| Peptide | `VVVGADGVGK` — KRAS G12D, residues 7–16 |
| MHC | HLA-A\*11:01 heavy chain (281 aa) + β2-microglobulin (99 aa) |
| Wild type | `VVVGAGGVGK` — differs at **position 6 only** |

The design problem is not "bind this surface". It is to discriminate a single
Gly→Asp substitution in a 10-mer lying in an MHC groove, where almost all of the
available surface is identical between target and off-target. Mutation
specificity therefore outranks affinity throughout, and the ranking is
lexicographic: **specificity → peptide-dependence → affinity**. A tight binder
that also binds wild type is worthless here; a weak one that discriminates is not.

Both sequences were verified rather than trusted. The heavy chain's first 180
residues match chain A of the design target exactly (180/180), which confirms the
allele and, more usefully, that the residue numbering used throughout this work is
the numbering the challenge specifies. β2m matches residues 1–99 of 7OW3 chain B;
the crystal carries an additional initiator Met at residue 0, so the challenge
gives the mature chain. Both checks are tests (`tests/test_constants.py`).

## 2. Design target construction (stages 03–04)

Structures surveyed: 7OW3, 7OW4, 7OW5, 7OW6, 7PB2.

7OW3 is the wild-type complex and **carries Gly at peptide position 6**, so Asp6
does not exist in any crystal structure of this complex at the conformer used
here. It was therefore built:

1. **Rotamer scan** (stage 02) over χ1 × χ2 for an Asp at position 6, scoring
   clashes, closest approach to the groove, and side-chain SASA
   (`results/rotamer_scan/`).
2. **Target construction** (stage 03) at **χ1 = −70°, χ2 = +120°** — a
   solvent-exposed rotamer, consistent with 7OW6 where a TCR contacts Asp6
   directly.
3. **Protonation and minimisation** (stage 04), then cropping to heavy-chain
   residues **1–180** (α1/α2) plus the peptide. α3 and β2m are removed because
   the binder never sees them.

Final target: `data/targets_relaxed/7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb`,
chain A = groove, chain C = peptide. An alternative rotamer
(χ1 = −180°, χ2 = 0°) was built and relaxed identically, and is used in §7 to
quantify how much the results depend on this choice.

## 3. Interface definition — Levy (2010) rSASA

Used **everywhere** interface residues are needed. No distance cutoff is used
anywhere in this work to define an interface; closest-approach distances are
reported descriptively and decide nothing.

Per `kras_pmhc.interface`:

1. Per-residue SASA computed twice — in the complex, and with partner atoms
   removed at identical coordinates (Shrake–Rupley, 500 points, 1.4 Å probe).
2. rSASA = SASA ÷ the **Tien et al. (2013) theoretical** maximum ASA. The
   theoretical rather than empirical maxima: the empirical values are drawn from
   observed structures and so already reflect burial, which is the quantity
   being measured.
3. ΔrSASA = rSASA_monomer − rSASA_complex, with a residue counted as changed on
   binding only if its absolute SASA loss exceeds 0.1 Å².
4. Classification at rSASA 0.25 → interior / surface / support / rim / core.
   **Interface = support + rim + core; buried interface = core + support.**

Applied over two *groups* of chains rather than chain-by-chain. On a three-chain
system the literal per-chain reading counts the HLA↔peptide interface, which
swamps the binder signal; single-chain groups reproduce the per-chain case
exactly.

**Validation against a known answer.** On 7OW6: peptide vs HLA puts Asp6 in
`core` (rSASA 0.572 → 0.089, 93.1 Å² lost), agreeing with an independent
side-chain SASA measurement of 7.4 Å²; peptide vs TCR puts the same residue in
`rim` (0.572 → 0.488, 16.2 Å² lost). The TCR therefore reads the aspartate
without desolvating it — recovering Poole et al.'s indirect-electrostatic
mechanism by a method that uses no distances at all. Both are asserted in
`tests/test_interface.py`.

## 4. Backbone generation (stage 05)

**RFdiffusion** via NVIDIA's hosted endpoint
(`biology/ipd/rfdiffusion/generate`), 50 designs, ~16 s each.

| Parameter | Value |
|---|---|
| `contigs` | `A1-180/0 C1-10/0 60-75` |
| `hotspot_res` | `C5, C6, C7` |
| `diffusion_steps` | 50 |
| noise scale | **endpoint default — not settable** |

Hotspots are centred on the mutation; C5/C6/C7 are Ala/**Asp**/Gly, so Asp6 is
the only side chain of substance in the set. Three hotspots is within the 3–6 the
authors recommend.

**Deviation, with consequences in §7:** the authors recommend `noise_scale = 0.5`
for protein–protein interfaces *because it improves designability*. The hosted
endpoint exposes no noise field. GPU routes that would have allowed it were
attempted first and failed for environmental reasons (documented in
`docs/colab_runbook.md`); this was the route that worked.

Three properties of the output had to be established empirically, because
assuming them produced plausible-looking but meaningless numbers:

* **The binder comes first**, not last as the contig string implies. Resolved per
  design by superimposing the putative target on the reference: 0.09 Å for
  binder-first against 21.70 Å for target-first.
* **Output is backbone-only** (N, CA, C, O) for the target as well as the binder,
  so Asp6 has no carboxylate. Since the target is held rigid, its side chains are
  known exactly — the stage-04 reference — and are reattached by splicing the
  generated chain onto it.
* **Output is recentred**, by ~81 Å RMS. Undone by superimposing on the target
  backbone. The shape check uses N/CA/C only, because the carbonyl O is rebuilt
  from idealised geometry (0.37 Å rms, 2.4 Å max, against ~0.1 Å for N/CA/C).

### Design gate

All five must hold. 20 of 50 designs pass criteria 1–4; **9 pass all five.**

| # | Criterion | Threshold | Basis |
|---|---|---|---|
| 1 | No backbone clash | ≥ 2.0 Å to target backbone | below any real heavy-atom contact |
| 2 | Asp6 is an interface residue | support/rim/core | Levy definition |
| 3 | Asp6 buried surface | **≥ 10 Å²** | gap in the data: …15.7, 14.7, 11.2, then 6.7, 4.3, 0 |
| 4 | Peptide share of target-side interface | **≥ 10%** | weak veto on HLA-dominated interfaces |
| 5 | Binder compactness | **Rg ≤ 1.7×** globular | designs split 1.13–1.62 vs 2.51–2.77, nothing between |

Criterion 1 is judged on the **backbone**, not all heavy atoms: RFdiffusion
conditions on backbone frames and cannot see the target's side-chain rotamers —
which this work chose — so it is not judged for hitting them. An all-atom rule
rejected 5 of 50 designs whose backbones were 2.5–5.9 Å clear and whose only
close contacts were flexible arginines. Those designs are flagged
`needs_repack` instead.

Criteria 3 and 4 encode the two leading objectives *separately and in order*. An
earlier version used criterion 4 alone at ≥ 20%, which inverted the project's own
ranking: it rejected a design burying 41 Å² of Asp6 (17% share) while admitting
one burying 15 Å² (37%) and one burying 6.7 Å² (41%) — the last being a
pan-peptide binder rather than a G12D reader.

Criterion 5 was added after generation, having been missing. 11 of the 20 designs
passing 1–4 were **single extended α-helices**, 88–106 Å end to end, which
satisfy every interface test — a helix laid across the groove buries plenty of
surface — while being unfoldable as standalone proteins. ProteinMPNN confirmed it
independently by answering them with poly-alanine.

## 5. Sequence design (stage 06)

**ProteinMPNN** via `biology/ipd/proteinmpnn/predict`. 10 sequences per surviving
backbone at `sampling_temp = 0.1`, designing the binder chain with the HLA and
peptide held as **fixed context** (the response confirms
`fixed_chains=['A','C'], designed_chains=['B']`), so the model is conditioned on
the target's sequence as well as its shape.

Two parameters were measured rather than assumed:

* `use_soluble_model = False`. The soluble model returned 71–80% hydrophobic
  sequences with runs of 43–53 alanines, against 43–63% for the default.
* `random_seed` **is not honoured** — two identical requests returned different
  sequences. Stage 06 is therefore not reproducible by re-running, so the raw
  multi-FASTA responses are preserved verbatim (`results/sequences/`) and the
  script will not overwrite one without `--force`.

90 sequences. Two exclusions leave **45 across 7 of 9 backbones**:

* **low complexity** — any residue type ≥ 40% (35 sequences). Not noise: it is
  what ProteinMPNN returns when a backbone has no core to pack.
* **a lone cysteine** (10 sequences, all with exactly one, so none could form a
  disulfide). An unpaired cysteine is a liability in a protein destined for
  expression.

Two backbones yield nothing usable. One of them, `nim3_0`, has the **best MPNN
score of the entire set** (0.683) and is 87% a single residue — score alone would
have selected the worst design.

## 6. Counter-screen (stages 07, 10, 08)

**Boltz-2** cofolding via NVIDIA's hosted endpoint, every selected sequence
against the G12D and wild-type peptides — which differ at position 6 alone, so
the pair isolates the mutation. 87 of 90 predictions completed; the run stalled
on a hung request, and the HTTP client sets no socket timeout.

Each prediction is **protonated and minimised** before scoring (stage 10).
PRODIGY counts atomic contacts across the interface, so its answer depends on
side-chain placement, and raw predicted side chains are the weakest part of that
input — the minimisations moved the potential energy from about −5,000 to
−14,400 kcal/mol, which is a measure of how far from a force-field minimum they
started. The protocol is stage 04's, and for stage 04's reason: amber14 with GBn2
implicit solvent, backbone heavy atoms restrained at 10 kcal/mol/Å², every side
chain free, no dynamics. The predicted pose is what is under test and must not
drift — observed backbone drift is ~0.09 Å — while the side chains, which are
what PRODIGY counts, are free to settle. **Both conditions of a given sequence
receive a byte-identical protocol**, which is the entire basis of the ΔΔG
comparison: relaxing them differently would make them differ for reasons
unrelated to position 6.

Cofolding is what makes PRODIGY possible at all: stage 05 produces backbones
only, and nothing contact-based can run on those.

**One result of this screen undermines its own premise, and is reported rather
than buried.**

In the cofolded models Asp6's rSASA within the pMHC is **0.11–0.12**, against
**0.59** in the crystal-derived target. The models bury the residue under test
inside the groove, where 7OW6 shows a TCR contacting it directly. Six of seven
designs consequently show zero peptide contact *by construction*. Likely causes,
in order: single-sequence mode (no MSA) degrading the peptide register, Boltz
modelling the decamer differently, and this work's own choice of an exposed
rotamer in §2.

What the screen does say, and it is not encouraging: the same predictions that
place the real HLA–peptide interface at **ipTM 0.94–0.96** place the binder at
**0.16–0.54** against the peptide and 0.15–0.59 against the HLA. Note that the
*overall* ipTM of 0.72–0.80 is carried almost entirely by the HLA–peptide
interface and would badly overstate the binders if quoted as a binder score.

G12D-minus-WT differences are ±0.07 with three of seven the wrong way round —
what one diffusion sample of noise produces.

A second chain-labelling trap was caught here: Boltz-2 **does not return the
chain ids it is given**, relabelling them sequentially by submission order.
Chains are resolved from their sequences, and the script refuses to proceed if
any cannot be identified.

## 7. Ranking and selection (stage 09)

Lexicographic, in the project's stated order of precedence:

1. **Mutation specificity** — **PRODIGY ΔΔG = ΔG(G12D) − ΔG(wild type)**, both
   from the same binder cofolded and minimised against peptides differing at
   position 6 alone. Negative favours the mutant. This is the primary criterion
   and what the submitted CSV is ordered by.
2. **Peptide-dependent binding** — the peptide's share of the target-side
   interface on the crystal-derived target, ≥ 10% for every design here.
3. **Affinity** — ProteinMPNN score, as the final tie-break.

Ranking is per **sequence**, not per backbone: each sequence has its own pair of
predictions, so each has its own ΔΔG.

**20 sequences are submitted**, the maximum the competition accepts from this
track. `molecule_class` is `protein` for all of them — single-chain de novo
miniproteins of 61–75 residues, which falls in the competition's *protein
minibinder* evaluation category (40–100 aa).

**Declared epitope.** Every design was generated with RFdiffusion hotspots on
peptide positions 5, 6 and 7, and selected for burying peptide residue 6. The
intended epitope is therefore the peptide and the groove immediately around it,
centred on position 6 — the one residue that differs between the G12D and
wild-type targets. This is stated so it need not be inferred from a prediction.

Two screens accompany each design in the CSV: **mean Kyte-Doolittle hydropathy**,
computed locally, and **complex pLDDT** from the cofolding response. The third
screen named by the competition, **ipSAE, could not be computed** — see §8.

`--rank-by burial` reproduces an alternative ordering on Asp6 buried area
instead of ΔΔG. The two disagree substantially, and §8 item 7 gives the reason
to treat the ΔΔG ordering with caution.

## 8. Stated shortcomings

Listed because they bear on how the results should be read, not as hedging.

1. **No design is validated as a binder, computationally or experimentally.**
   The one orthogonal check that ran placed binder interface confidence at
   0.16–0.54 where it placed a real interface at 0.94–0.96.
2. **The counter-screen could not test the intended interaction** (§6). There is
   therefore **no evidence of mutation specificity** — only of geometric
   engagement with a modelled Asp6.
3. **Asp6 is constructed, not observed.** Re-scoring all 50 designs against an
   independently built and relaxed alternative rotamer moves **9 of 50 gate
   verdicts** (3 lost, 6 gained; 17 pass under both). Individual swings are
   large: one design falls from 36.8 Å² of Asp6 burial to 4.0 Å². Selection uses
   the canonical rotamer by decision; `rotamer_robust` in `designs.tsv` records
   which designs survive the alternative, and `results/rotamer_sensitivity.tsv`
   has the full table.
4. **All backbones were sampled at default noise** (§4), which the authors'
   own guidance says costs designability — and the extended-helix problem in §4
   is what that cost looks like.
5. **Every binder is helical**, 84–97% with zero strand content across all 50
   designs. The topology search never explored sheets.
6. **No sequence was validated by refolding.** An AF2/Boltz self-consistency
   check of designed sequence against intended backbone was planned and not run.
7. **The ranking criterion is, by its own evidence, noise-dominated.** All
   sequences on one backbone share the same intended pose, so a genuine Asp6
   specificity signal should be similar across them. It is not. The spread of
   ΔΔG *within* a backbone is **sd 1.23 kcal/mol** while the spread *between*
   backbone means is **sd 0.71** — a signal-to-noise ratio of 0.58. Across all
   43 sequences with both conditions, 18 favour G12D where a coin flip predicts
   21.5. Differences of ~1 kcal/mol between designs are therefore inside the
   noise, and selecting the most negative values selects partly for noise. The
   submission is ordered by ΔΔG because that is the project's stated criterion
   and the best available; this item is the honest measure of how much weight the
   order can carry.
8. **No apo condition was predicted.** The competition's second-hardest stated
   goal is no binding to the HLA with no peptide, and this submission carries no
   evidence on it. The condition was staged for a Colab run that was not
   performed, and the hosted endpoint stopped serving requests before it could
   be done there.
9. **The cofolded target is a truncation.** Heavy-chain residues 1–180 (α1/α2)
   plus the peptide, without β2-microglobulin or α3, where the competition
   specifies modelling the whole complex. The groove face the binder was
   designed against is present, but the ΔG values are not directly comparable
   with numbers computed on the full construct.
10. **Only two of the three named screens are reported.** ipSAE could not be
    computed: the hosted endpoint returns the PAE matrix as `null`, and ipSAE is
    derived from PAE. A local Boltz run would supply it.
11. **Three of 90 cofolding predictions are missing**, so two sequences lack a
    ΔΔG and sort last.
12. **Five designs carry `needs_repack`** — a target side chain sits within 2.0 Å
    of binder backbone with the backbone itself clear. Minimisation will have
    relieved some of that, but it was not verified per design.
13. **Scale.** 50 backbones is small for de novo binder design, where hit rates
    are low and campaigns routinely sample thousands.

## 9. Reproducibility

```bash
conda activate kras-pmhc
make fetch survey rotamers targets relax   # stages 00-04
make nim-key                               # store the API key (hidden prompt)
make nim                                   # stage 05: 50 backbones
make backbones                             # scoring and the five-criterion gate
make rotamer-sensitivity                   # how far the gate depends on the rotamer
make sequences select                      # stage 06: ProteinMPNN, then filtering
make cofold PER_DESIGN=10                  # stage 07: cofold against G12D and WT
make minimise-cofold                       # stage 10: protonate + minimise
make score-cofold-min                      # stage 08: Levy + ipTM + PRODIGY
make package                               # stage 09 -> submission/
make test lint
```

Every threshold quoted here is a named constant in the script that applies it,
and every number is traceable to a file:

| | |
|---|---|
| `data/backbones/MANIFEST.tsv` | all 50 designs, every gate input and verdict |
| `results/interface/` | per-residue Levy tables behind each interface number |
| `results/rotamer_sensitivity.tsv` | the same designs scored on the alternative rotamer |
| `results/sequences/` | raw ProteinMPNN responses, verbatim |
| `data/cofolding_nim/` | predictions, minimised models, and minimisation energies |
| `results/cofolding/SCORES_min.tsv` | the scores the ranking is computed from |
| `submission/` | the CSV, FASTA, evidence table, this document, and structures |

The repository has 153 tests. The analysis code is covered directly; the two GPU
notebooks are checked for syntax and for pinned parameters, because a typo in one
only surfaces after minutes of installation on a session that may not be granted
again.

**Stages that were not fully scripted.** Backbone generation and cofolding both
needed a GPU this machine does not have. Both ended up on hosted endpoints and
are scripted; the Colab notebooks written for them, and the manual steps they
would have required, are recorded in `docs/colab_runbook.md` along with the
environment failures that made them necessary.

## References

- Levy, E.D. (2010) *J Mol Biol* **403**:660 — rSASA interface classification.
- Tien, M.Z. et al. (2013) *PLoS ONE* **8**:e80635 — theoretical maximum ASA.
- Watson, J.L. et al. (2023) *Nature* **620**:1089 — RFdiffusion.
- Dauparas, J. et al. (2022) *Science* **378**:49 — ProteinMPNN.
- Wohlwend, J. et al. (2024/2025) — Boltz-1 / Boltz-2.
- Xue, L.C. et al. (2016) *Bioinformatics* **32**:3676 — PRODIGY.
- Kyte, J. & Doolittle, R.F. (1982) *J Mol Biol* **157**:105 — hydropathy scale.
- Eastman, P. et al. (2017) *PLoS Comput Biol* **13**:e1005659 — OpenMM.
- Poole, A. et al. (2022) *Nat Commun* **13**:5333 — TCR recognition of KRAS
  G12D on HLA-A\*11:01.
