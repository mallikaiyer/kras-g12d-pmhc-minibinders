<!-- Kept in the repository, deliberately not part of the submission. -->
<!-- The submitted methods paper is docs/methods.md. -->

# Stated shortcomings

Listed because they bear on how the results should be read.

1. **No design is validated as a binder, computationally or experimentally.** The
   one orthogonal check that ran placed binder interface confidence at 0.16–0.54
   where the same predictions placed a real interface at 0.94–0.96.
2. **The counter-screen does not present the residue under test.** In the
   cofolded models, peptide residue 6 has a relative SASA within the pMHC of
   0.11–0.12, against 0.59 in the prepared crystal-derived target — the
   predictions bury the aspartate inside the groove, where crystallography shows
   a receptor contacting it directly. Six of seven designs consequently show zero
   peptide contact *by construction*. Likely causes, in order: single-sequence
   prediction degrading the peptide register, the predictor modelling a decamer
   differently, and this work's own choice of a solvent-exposed rotamer. The
   screen is therefore **uninformative about the intended interaction rather than
   evidence against it** — and, correspondingly, there is **no positive evidence
   of mutation specificity**, only of geometric engagement with a modelled
   aspartate.
3. **The modelled aspartate is not observed.** The free complex used as the
   design template carries wild-type glycine at position 6. Nine of 50 gate
   verdicts move on changing the rotamer (§3.6).
4. **The open-state population is not measured.** The approach rests on
   conformational selection, which requires the open state to be populated enough
   in the free ensemble. No molecular dynamics or free-energy calculation was
   run; the structural evidence establishes only that the open conformer is
   crystallographically observable. If the open state is rare, the −RT ln *p*
   penalty is larger than assumed.
5. **"The open conformation" is an ensemble average treated as a state**, and the
   backbone used is from the wild-type complex — the free G12D peptide is, if
   anything, more mobile.
6. **No apo condition was predicted.** Objective 2 — no binding to the
   peptide-free complex — carries no evidence in this submission.
7. **The cofolded target is a truncation**: heavy-chain residues 1–180 plus
   peptide, without β2-microglobulin or α3. The groove face the designs were
   built against is present, but the free energies are not directly comparable
   with values computed on the full construct.
8. **All backbones were sampled at the generator's default noise**, which the
   method's own guidance says costs designability; the extended-helix problem in
   §3.6 is what that cost looks like.
9. **Every binder is helical** — 84–97% helix with zero strand content across
   all 50 designs. The topology search never explored sheets.
10. **No sequence was validated by refolding.** A self-consistency check of
    designed sequence against intended backbone was planned and not run.
11. **One diffusion sample per prediction**, so no error bars on any
    cofolding-derived quantity.
12. **Only two of the three named screens are reported** (§3.10).
13. **Three of 90 predictions are missing**, so two sequences lack a ΔΔG and sort
    last.
14. **Five designs had a target side chain within 2.0 Å of binder backbone**
    before minimisation, with the backbone itself clear. Minimisation will have
    relieved some of that; it was not re-verified per design.
15. **Scale.** 50 backbones is small for de novo binder design, where hit rates
    are low and campaigns routinely sample thousands.

These bear on how the results should be read. They are recorded here
rather than in the submitted paper because the submission is a ranked CSV
and a methods description; a self-assessment of the work's weaknesses is
for the record, not for the submission.

# Potential improvements

Each addresses one or more of the shortcomings above, roughly in order of how
much it would change the result.

1. **Attempt to fix the counter-screen so it presents the residue under test**
   (2, 11). The single highest-value change. Options, cheapest first: cofold
   with a multiple-sequence alignment rather than in single-sequence mode; take
   several diffusion samples and check the peptide register in each; constrain
   or template the peptide to the open backbone.
2. **Measure the open-state population** (4, 5). The conformational-selection
   premise is currently supported by crystallographic observability alone.
   Enhanced-sampling molecular dynamics on the free G12D complex would put a
   number on *p*, and therefore on the penalty the approach is paying. If the
   open state turns out to be rare, the premise needs revisiting rather than
   refining.
3. **Add the apo condition** (6). Objective 2 carries no evidence at all.
   Cofolding each design against the peptide-free groove is the same machinery
   already built, one more condition — it was dropped for time.
4. **Regain control of the generator's noise scale, and scale up** (8, 15).
   Running RFdiffusion via Colab or Kaggle rather than through the NVIDIA hosted
   endpoint restores `noise_scale_ca`/`noise_scale_frame` = 0.5, which the
   authors recommend for interfaces precisely because it improves designability,
   and removes the per-call ceiling that capped this run at 50 backbones.
   Hundreds to thousands is the published norm.
5. **Add a refolding self-consistency filter** (10). Predict each designed
   sequence alone and require the result to match the backbone it was designed
   for. This catches sequences that score well without encoding their intended
   fold, which nothing here tested.
6. **Explore non-helical topologies** (9). Every design is 84–97% helix with
   zero strand content. Secondary-structure or fold conditioning, or the beta
   checkpoint, would at least test whether a sheet-containing binder reads this
   epitope better.
7. **Carry a rotamer ensemble rather than one rotamer** (3). Nine of 50 gate
   verdicts move between two rotamers. Scoring against several and requiring a
   design to pass under all of them converts a modelling assumption into a
   robustness criterion. The data for two rotamers already exists in
   `results/rotamer_sensitivity.tsv`.
8. **Score on the full construct** (7). Including β2-microglobulin and α3 makes
   the predicted free energies comparable with literature values computed on the
   intact complex, at the cost of a larger prediction.
9. **Report ipSAE** (12). The hosted endpoint returns a null PAE matrix and
   ipSAE is derived from it; a local prediction, or an endpoint that returns PAE,
   makes the third named screen computable.
10. **Re-verify clashes after minimisation** (14). Cheap, and closes a gap that
    is currently stated rather than resolved.
