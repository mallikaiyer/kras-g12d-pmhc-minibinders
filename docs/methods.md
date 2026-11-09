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

## 6. Counter-screen (stages 07–08)

**Boltz-2** cofolding, each selected sequence against the G12D and wild-type
peptides — which differ at position 6 alone — to isolate the mutation. PRODIGY
on the resulting complexes for binding free energy. Cofolding is what makes
PRODIGY possible at all: it counts atomic contacts, so it needs side chains on
both partners, and stage 05 produces backbones only.

**This screen did not work, and its output is reported but not ranked on.**

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

Lexicographic, on the **crystal-derived** target rather than the cofolded models,
for the reason in §6 — ranking on structures that do not present the residue
under test would rank on an artefact.

1. **Mutation specificity** — Asp6 buried surface area (≥ 10 Å² for all).
2. **Peptide-dependent binding** — peptide share of target-side interface
   (≥ 10% for all).
3. **Affinity** — ProteinMPNN score, used only to choose between sequences on the
   same backbone.

10 sequences are submitted, allocated **round-robin across all 7 backbones**
before any backbone gets a second slot. Submitting several sequences from one
backbone is a bet on that backbone, and nothing here justifies concentrating the
set.

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
7. **One diffusion sample per cofolding prediction**, so no error bars; the ±0.07
   specificity differences in §6 are not distinguishable from noise.
8. **Five designs carry `needs_repack`** — a target side chain sits within 2.0 Å
   of binder backbone with the backbone itself clear. Any energetic claim for
   those requires repacking first.
9. **Scale.** 50 backbones is small for de novo binder design, where hit rates
   are low and campaigns routinely sample thousands.

## 9. Reproducibility

```bash
conda activate kras-pmhc
make fetch survey rotamers targets relax   # stages 00-04
make nim                                   # stage 05 generation
make backbones                             # scoring and the design gate
make sequences select                      # stage 06
python scripts/07b_nim_cofolding.py        # stage 07 counter-screen
python scripts/08_score_cofolded.py        # stage 08
python scripts/09_rank_and_package.py      # stage 09 -> submission/
make test lint
```

Every threshold quoted here is a named constant in the script that applies it,
and every number is traceable to a file: per-residue Levy tables in
`results/interface/`, raw service responses in `results/sequences/` and
`data/cofolding_nim/`, and the full design manifest in
`data/backbones/MANIFEST.tsv`.

## References

- Levy, E.D. (2010) *J Mol Biol* **403**:660 — rSASA interface classification.
- Tien, M.Z. et al. (2013) *PLoS ONE* **8**:e80635 — theoretical maximum ASA.
- Watson, J.L. et al. (2023) *Nature* **620**:1089 — RFdiffusion.
- Dauparas, J. et al. (2022) *Science* **378**:49 — ProteinMPNN.
- Wohlwend, J. et al. (2024/2025) — Boltz-1 / Boltz-2.
- Xue, L.C. et al. (2016) *Bioinformatics* **32**:3676 — PRODIGY.
- Poole, A. et al. (2022) *Nat Commun* **13**:5333 — TCR recognition of KRAS
  G12D on HLA-A\*11:01.
