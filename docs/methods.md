# De novo minibinders against the KRAS G12D neoepitope on HLA-A\*11:01

Anthropic × Adaptyv 2026 protein design competition, Challenge 5 (peptide–MHC).
Solo entry, Track 3.

Accompanies a ranked CSV of 20 designed sequences.

---

## 1. The challenge

The target is the KRAS G12D neoepitope `VVVGADGVGK` — residues 7–16 of mutant
KRAS — presented on HLA-A\*11:01 with β2-microglobulin. The off-targets are the
wild-type peptide `VVVGAGGVGK` on the same HLA, and the same HLA with no peptide
at all.

The objectives, hardest first:

1. **No binding to the wild-type peptide complex.**
2. **No binding to the peptide-free complex.**
3. **Strong binding to the G12D complex.**

The difficulty is entirely in the first objective. Target and off-target differ
by a single methylene-to-carboxylate substitution at one position of a ten-residue
peptide, and the peptide itself contributes a minority of the solvent-accessible
surface of the complex. Nearly every atom a binder could contact is identical
between the two.

Affinity is the *least* important objective, and over-optimising it is actively
counterproductive. The therapeutic precedent makes this concrete: the parent JDI
TCR against this epitope binds G12D at 63 µM with no detectable wild-type
binding, and as it was affinity-matured through 25.8 µM and 3.0 µM to 425 nM,
**wild-type binding appeared**. Micromolar affinity with clean discrimination
beats nanomolar affinity without it.

## 2. Approach and premise

### The discriminating feature

Peptide position 6 is the only chemical difference between the target and the
wild-type off-target. The aspartate presents a **carboxylate**: a charged group
a designed binder can reach with a salt bridge, or satisfy with backbone donors.
Wild-type glycine presents no side chain at all.

The discriminator available to a binder is therefore **aspartate versus
nothing** — not two competing chemistries to be told apart, but presence against
absence.

This is unusually favourable for negative design. A pocket built to present a
cation and backbone donors sized for a carboxylate is, against the wild-type
peptide, a desolvated and unsatisfied cation — *destabilising*, not merely
unrewarded. The same structural feature serves objectives 1 and 3 simultaneously.
Objective 2 follows for free: the peptide-free groove has nothing for such a
pocket to grab.

### Conformational selection on the open state

Crystallography of this system shows the peptide centre adopting two
arrangements. In the free complex the mid-peptide backbone is **open**, with
position 6 projecting out of the groove; in the receptor-bound complex it is
**closed**, with position 6 turned down into the groove and engaged by the HLA
F-pocket.

This work designs against the **open** backbone, and the mechanistic claim is
**conformational selection rather than induced fit**: the open state need only be
*populated* in the free ensemble, not dominant. A binder that engages it shifts
the equilibrium toward it. The cost is a penalty on apparent affinity of
−RT ln *p* — about 1.4 kcal/mol if the open state is 10% populated, 2.7 kcal/mol
at 1%.

Four observations support the choice:

1. **The disorder in the free G12D structure is evidence for the premise, not
   against it.** No density is modelled at peptide residues 5–6 in any of the
   four copies of the free G12D complex. Conformational selection requires a
   populated alternative state, and disorder is what that looks like in a
   diffraction experiment.
2. **The open state is an observation, not a model.** Two independent copies of
   the free wild-type complex agree to within 0.5 Å, and a wild-type
   receptor-bound control shows the open-to-closed flip is driven by receptor
   binding rather than by the mutation.
3. **Reachability was measured rather than assumed.** On that backbone, every
   clash-free aspartate rotamer leaves the side chain solvent-exposed at
   ≥ 38 Å².
4. **Negative design is built in**, as above.

### Known headwinds

Stated at the outset, because they bound what this approach can achieve.

The groove rewards burial: only the G12D peptide forms a stabilising network
with the HLA F-pocket in the closed state, so the closed conformation is *more*
stabilised for the mutant than for wild type — the design is working against the
mutant's own preference. And the entropic cost is paid twice: ordering a mobile
region, in a non-preferred conformation. For scale, the natural TCR paid
TΔS° = −8.13 kcal/mol across ~990 Å² of buried surface; a 60–75 residue
miniprotein buries roughly 800–900 Å², so there is margin but not much.

The expected outcome of this approach is a binder with a defensible mechanism and
a genuine chance at objective 1 — not a high-probability binder.

## 3. Methods

### 3.1 Structural basis

Five crystal structures of this epitope system were surveyed: the free wild-type
and free G12D complexes, two receptor-bound complexes (wild type and G12D with
the same engineered TCR), and a third bound complex with the parent TCR. The free
wild-type complex provided the open backbone used throughout; its peptide
position 6 is wild-type glycine, so the aspartate had to be modelled.

### 3.2 Aspartate rotamer construction

Position 6 was mutated glycine→aspartate and its side chain built by a scan over
χ1 × χ2, scoring each rotamer for steric clash, closest approach to the groove,
and side-chain solvent accessibility. **χ1 = −70°, χ2 = +120°** was selected: a
clash-free, solvent-exposed rotamer consistent with the receptor-bound structure
in which the aspartate is contacted directly.

An independent alternative, **χ1 = −180°, χ2 = 0°**, was built and carried
through identically for the sensitivity analysis in §3.6.

### 3.3 Target preparation

Hydrogens were added at pH 7.4 and the structure energy-minimised with **amber14**
and **GBn2 implicit solvent**, backbone heavy atoms restrained harmonically at
**10 kcal/mol/Å²** and all side chains free. No dynamics were run. The restraint
preserves the experimentally observed backbone — the central claim of the
approach — while letting the environment respond to the new aspartate.

The design target was then cropped to **heavy-chain residues 1–180** (the α1/α2
groove) **plus the peptide**. α3 and β2-microglobulin were removed: a binder
engaging the peptide and the groove rim does not contact them.

### 3.4 Interface definition

Every interface in this work is identified by the **Levy (2010) relative-SASA
classification**. No distance cutoff is used anywhere to define an interface;
closest-approach distances are reported descriptively and decide nothing.

Per-residue SASA is computed twice — in the complex, and with partner atoms
removed at identical coordinates (Shrake–Rupley, 500 points, 1.4 Å probe — then
divided by the **Tien et al. (2013) theoretical** maximum accessible area. The
theoretical rather than empirical maxima: empirical values are derived from
observed structures and so already embed burial, which is the quantity being
measured. A residue counts as changed on binding only if its absolute SASA loss
exceeds 0.1 Å². Classification at relative SASA 0.25 gives interior, surface,
support, rim and core, with **interface = support + rim + core**.

The classification is applied over two *groups* of chains rather than
chain-by-chain, because on a three-chain system the per-chain reading counts the
HLA–peptide interface and swamps the binder signal.

**Validation against a known answer.** On the G12D receptor-bound structure, the
peptide-versus-HLA grouping places Asp6 in **core** (relative SASA 0.572 → 0.089,
93.1 Å² lost), agreeing with an independent side-chain SASA measurement. The
peptide-versus-TCR grouping places the same residue in **rim** (0.572 → 0.488,
16.2 Å² lost): the receptor reads the aspartate without desolvating it. That
recovers the published indirect-electrostatic mechanism for this system by a
method that uses no distances at all.

### 3.5 Backbone generation

**RFdiffusion**, via a hosted inference endpoint, 50 designs.

| Parameter | Value |
|---|---|
| Contigs | `A1-180/0 C1-10/0 60-75` (target held fixed, 60–75 residue binder) |
| Hotspots | peptide positions **5, 6, 7** |
| Diffusion steps | 50 |
| Noise scale | endpoint default — **not settable** |

Hotspots are centred on the mutation. Positions 5/6/7 are Ala/**Asp**/Gly, so the
aspartate is the only side chain of substance in the set, and three hotspots sits
within the three to six the method's authors recommend.

The noise scale is a known deviation, and its consequences are visible in
§3.6: the authors
recommend 0.5 for protein–protein interfaces *because it improves designability*,
and the endpoint used exposes no such field, so all 50 designs were sampled at
the default.

Three properties of the generator output were established empirically, because
assuming them produced plausible-looking but meaningless interface numbers:

- **The generated chain is written first**, not last as the contig string
  implies. Resolved per design by superimposing the putative target onto the
  reference: 0.09 Å under the correct assignment against 21.70 Å under the
  wrong one.
- **Output carries backbone atoms only** (N, Cα, C, O) for the target as well as
  the binder, so the aspartate has no carboxylate. Because the target is held
  rigid during generation its side chains are known exactly, and were restored
  by splicing the generated chain onto the prepared target.
- **Output is recentred**, by ~81 Å RMS, undone by superposition on the target
  backbone. Shape agreement is assessed over N/Cα/C only: the carbonyl oxygen is
  rebuilt from idealised geometry and differs by 0.37 Å RMS and up to 2.4 Å,
  against ~0.1 Å for the frame atoms.

### 3.6 Design gate

Five criteria, all required. **9 of 50 designs pass.**

| # | Criterion | Threshold | Basis for the threshold |
|---|---|---|---|
| 1 | No backbone clash | ≥ 2.0 Å to target backbone | below any real heavy-atom contact |
| 2 | Asp6 is an interface residue | support, rim or core | Levy definition |
| 3 | Asp6 buried surface area | **≥ 10 Å²** | a gap in the observed data: …15.7, 14.7, 11.2, then 6.7, 4.3, 0 |
| 4 | Peptide share of target-side interface | **≥ 10%** | weak veto on HLA-dominated interfaces |
| 5 | Binder compactness | **Rg ≤ 1.7 ×** globular | designs split 1.13–1.62 against 2.51–2.77, nothing between |

Criterion 1 is judged on the target **backbone**, not all heavy atoms, because
the generator conditions on backbone frames and cannot see the target's
side-chain rotamers — which this work chose. An all-atom rule rejected five
designs whose backbones were 2.5–5.9 Å clear and whose only close contacts were
flexible arginines.

Criteria 3 and 4 encode the two leading objectives separately and in order. An
earlier form of the gate used criterion 4 alone at ≥ 20% and inverted that
priority: it rejected a design burying 41 Å² of Asp6 at 17% share while admitting
one burying 15 Å² at 37%, and one burying 6.7 Å² at 41% — the last being a
pan-peptide binder rather than a neoepitope reader.

Criterion 5 was added after generation, having been absent. Eleven of the twenty
designs that passed criteria 1–4 were **single extended α-helices**, 88–106 Å end
to end, which satisfy every interface test — a helix laid across the groove
buries plenty of surface — while being unfoldable as standalone proteins.
Sequence design confirmed it independently by answering them with poly-alanine.

**Rotamer sensitivity.** Because the aspartate is modelled rather than observed,
all 50 designs were re-scored against the independently built alternative
rotamer, changing nothing else. **Nine of 50 gate verdicts move** (three lost,
six gained; 17 pass under both). Individual swings are large: one design falls
from 36.8 Å² of Asp6 burial to 4.0 Å². Selection uses the primary rotamer; the
submitted CSV records which designs survive the alternative.

### 3.7 Sequence design

**ProteinMPNN**, via a hosted endpoint. Ten sequences per surviving backbone at
sampling temperature **0.1**, designing the binder chain with the HLA and peptide
held as **fixed context**, so the model is conditioned on the target's sequence
as well as its shape.

Two parameters were set by measurement rather than default:

- The **soluble-trained model was rejected**: on a test backbone it returned
  sequences 71–80% hydrophobic with runs of 43–53 alanines, against 43–63% for
  the standard model.
- The **random seed is not honoured** by the endpoint — two identical requests
  returned different sequences — so this step is not reproducible by re-running,
  and the raw responses were retained as the record.

90 sequences were generated. Two filters were applied:

- **Low complexity**, any residue type at ≥ 40% of the sequence: 35 sequences.
  This is not sampling noise but what the method returns when a backbone has no
  core to pack.
- **A lone cysteine**: 10 sequences, each with exactly one, so none could form an
  intramolecular disulfide. An unpaired cysteine is a liability in a protein
  intended for expression.

**45 sequences across 7 backbones** survived. Two backbones yielded nothing
usable; one of them produced the best sequence score of the entire set (0.683)
while being 87% a single residue and 88% hydrophobic, so score alone would have
selected the worst design available.

### 3.8 Counter-screen by cofolding

**Boltz-2**, via a hosted endpoint, every surviving sequence against the G12D and
wild-type peptides — differing at position 6 alone, so the pair isolates the
mutation. 87 of 90 predictions completed. Single diffusion sample, 50 sampling
steps, 3 recycling steps, step scale 1.638, **single-sequence mode with no
multiple-sequence alignment**.

Cofolding is what makes contact-based scoring possible at all: generation
produces backbones without side chains, and nothing contact-based can run on
those.

Chain identity in the output is **resolved from the sequences**, not from the
returned chain labels: the endpoint relabels chains sequentially by submission
order rather than preserving the identifiers it is given, and trusting those
labels placed the binder where the peptide was expected.

### 3.9 Minimisation of the predicted complexes

Each cofolded complex was protonated and minimised before scoring, under the
protocol of §3.3 and for the same reason: amber14 with GBn2 implicit solvent,
backbone heavy atoms restrained at 10 kcal/mol/Å², all side chains free, no
dynamics.

Contact-based scoring depends on side-chain placement, and raw predicted side
chains are the weakest part of that input: minimisation moved the potential
energy from roughly −5,000 to −14,400 kcal/mol, which measures how far from a
force-field minimum they began. Backbone drift under the restraint was ~0.09 Å,
so the predicted pose — the thing under test — is preserved.

**Both conditions of a given sequence received an identical protocol.** That is
the entire basis of the differential in §3.11: relaxing them differently would
make them differ for reasons unrelated to position 6.

### 3.10 Scoring

On each minimised complex:

- **PRODIGY** predicted binding free energy, binder against the HLA-plus-peptide
  group, matching the Levy grouping used for every interface number here.
- **Levy classification** of the interface, including the class and buried area
  of peptide residue 6.
- **Interface confidence** from the cofolding response, read **per chain pair**.
  The overall interface-pTM of 0.72–0.80 is carried almost entirely by the real
  HLA–peptide interface and would badly overstate the binders if quoted as a
  binder score; the binder-to-peptide and binder-to-HLA pairs are the relevant
  figures.

Two of the three sequence and fold screens named by the competition accompany
each design: **mean Kyte–Doolittle hydropathy**, computed directly, and **complex
pLDDT** from the cofolding response. **ipSAE could not be computed**: the
endpoint returns the predicted aligned error matrix as null, and ipSAE is derived
from it.

### 3.11 Ranking and selection

Lexicographic, in the objectives' order of precedence:

1. **Mutation specificity** — **ΔΔG = ΔG(G12D) − ΔG(wild type)**, both from the
   same binder cofolded and minimised against peptides differing at position 6
   alone. Negative favours the mutant. This is the primary criterion and the
   order of the submitted CSV.
2. **Peptide-dependent binding** — the peptide's share of the target-side
   interface on the prepared crystal-derived target; ≥ 10% for every design
   submitted.
3. **Affinity** — sequence design score, as the final tie-break only.

Ranking is per **sequence** rather than per backbone, since each sequence has its
own pair of predictions and therefore its own ΔΔG.

**20 sequences are submitted**, the maximum this track accepts. All are
single-chain de novo miniproteins of 61–75 residues — `molecule_class` =
`protein`, falling in the competition's protein-minibinder category (40–100 aa).

**Declared epitope.** Every design was generated with hotspots on peptide
positions 5, 6 and 7 and selected for burying peptide residue 6. The intended
epitope is the peptide and the groove immediately around it, centred on position
6 — the one residue that differs between the G12D and wild-type targets. It is
stated here so that it need not be inferred.

## 4. References

- Levy, E.D. (2010) A simple definition of structural regions in proteins and its
  use in analyzing interface evolution. *J Mol Biol* **403**:660.
- Tien, M.Z., Meyer, A.G., Sydykova, D.K., Spielman, S.J. & Wilke, C.O. (2013)
  Maximum allowed solvent accessibilities of residues in proteins. *PLoS ONE*
  **8**:e80635.
- Shrake, A. & Rupley, J.A. (1973) Environment and exposure to solvent of protein
  atoms. Lysozyme and insulin. *J Mol Biol* **79**:351.
- Watson, J.L. et al. (2023) De novo design of protein structure and function
  with RFdiffusion. *Nature* **620**:1089.
- Dauparas, J. et al. (2022) Robust deep learning–based protein sequence design
  using ProteinMPNN. *Science* **378**:49.
- Wohlwend, J. et al. (2024) Boltz-1: democratizing biomolecular interaction
  modeling; and Passaro, S. et al. (2025) Boltz-2: towards accurate and efficient
  binding affinity prediction.
- Xue, L.C., Rodrigues, J.P., Kastritis, P.L., Bonvin, A.M. & Vangone, A. (2016)
  PRODIGY: a web server for predicting the binding affinity of protein–protein
  complexes. *Bioinformatics* **32**:3676.
- Kyte, J. & Doolittle, R.F. (1982) A simple method for displaying the hydropathic
  character of a protein. *J Mol Biol* **157**:105.
- Eastman, P. et al. (2017) OpenMM 7: rapid development of high performance
  algorithms for molecular dynamics. *PLoS Comput Biol* **13**:e1005659.
- Maier, J.A. et al. (2015) ff14SB: improving the accuracy of protein side chain
  and backbone parameters from ff99SB. *J Chem Theory Comput* **11**:3696.
- Nguyen, H., Roe, D.R. & Simmerling, C. (2013) Improved generalized Born solvent
  model parameters for protein simulations. *J Chem Theory Comput* **9**:2020.
- Poole, A. et al. (2022) Therapeutic high affinity T cell receptor targeting a
  KRAS G12D cancer neoantigen. *Nat Commun* **13**:5333.
