# Challenge 5 — KRAS G12D pMHC: design plan

**Target:** `VVVGADGVGK` (KRAS G12D, P01116 res 7–16) on HLA-A\*11:01 + β2-microglobulin
**Counter-targets:** wild-type `VVVGAGGVGK` on the same HLA; peptide-free HLA-A\*11:01/β2m
**Deadline:** Sun 8 Nov 2026, 23:59 AoE · **Rewritten 3 Nov 2026** · Track 3, solo

Ranking is lexicographic: (1) rejects wild-type, (2) rejects peptide-free HLA, (3) affinity. The winner must meet every objective; one honourable mention per objective. Affinity is the *last* tie-break, so a weak clean binder outranks a strong promiscuous one.

## Constraints

| Constraint | Consequence |
|---|---|
| MacBook Air M2, no CUDA | GPU work goes to free Colab T4 / AlphaFold Server. Local work is CPU only. |
| No Claude/Modal credits, no commercial licences | No Rosetta, no PyRosetta pipelines. PRODIGY is the affinity proxy. |
| ~2 h/day, 3–8 Nov | ~12 h of attention. Generation must run unattended. |
| Track 3, solo | ≤20 designs; only the first 20 count, in CSV row order; ≤10 per submission tested. **Rows 1–10 are the real submission.** |

**Fix first:** the mambaforge numpy is broken (`libcblas.3.dylib` not found). Every tool below depends on it.

---

## What the structures actually say

All five structures are from Poole et al., *Nat Commun* 13:5333 (2022), [10.1038/s41467-022-32811-1](https://doi.org/10.1038/s41467-022-32811-1).

| PDB | Contents | Residue 5–6 density |
|---|---|---|
| **7OW3** | WT pMHC, **free** (4 copies in ASU) | complete in peptide chains **C** and **L**; missing in F and I |
| **7OW4** | G12D pMHC, **free** (4 copies) | **absent in all four copies** |
| 7OW5 | WT pMHC + JDIa41b1 TCR | ordered |
| **7OW6** | G12D pMHC + JDIa41b1 TCR | ordered, real Asp6 |
| 7PB2 | G12D pMHC + parent JDI TCR | — |

**Position 6 changes conformation on receptor binding.** Measured distances to HLA Arg114 CZ (the residue Asp6 salt-bridges when buried):

| state | Cα → R114 | Cβ → R114 | side chain points |
|---|---|---|---|
| 7OW3 chain C (free) | 10.9 Å | 11.9 Å (virtual) | **away** from the groove |
| 7OW3 chain L (free) | 10.5 Å | 11.4 Å (virtual) | **away** from the groove |
| 7OW6 (TCR-bound) | 6.4 Å | 5.7 Å (real) | **into** the groove |

Two things change: the backbone drops ~4.3 Å, and the sign of (Cβ − Cα) inverts. In the **free** complex — which is what the assay presents — residue 6 projects away from the groove and is solvent-accessible. In the TCR-bound complex it is buried: the carboxylate sits 2.95 Å from HLA R114 and 3.08 Å from Q70, i.e. fully engaged by the HLA and unreachable by a binder.

**Backbone compatibility (checked).** φ/ψ at residue 6: 7OW3 C = −69.6°/+48.3°, 7OW3 L = −86.4°/+41.7°, 7OW6 = −37.2°/+117.7°. Both free copies have **negative φ**, so the open backbone is not Gly-exclusive (Gly-only territory is overwhelmingly positive φ) and an Asp can occupy it. The free copies sit in the bridge region between αR and β — allowed but not most-favoured, so relaxation is warranted. Note the *experimentally observed* Asp6 in 7OW6 is also in a sparsely populated φ/ψ combination: marginal dihedrals at position 6 are a property of this peptide in this groove, not an artefact of a model.

**The TCR's discrimination is not contact-mediated.** Min distance from the Asp6 carboxylate to the residues carrying the largest per-residue binding-energy differences in 7OW6:

| TCR residue | ΔΔG contribution | distance to D6 carboxylate |
|---|---|---|
| R50α | 2.0 kcal/mol | **13.5 Å** |
| K70α | 0.88 kcal/mol | **15.9 Å** |
| K94β | — | 13.2 Å |
| Q99β | — | 9.2 Å |
| H100β | — | 7.0 Å |
| N101β | — | 9.3 Å |

The paper's clause "directly positioned above the mutant residue" refers to **the CDR3β loop**, not to these residues individually; K70α is explicitly described as "further from the mutation site." The mechanism is long-range electrostatic complementarity to a more negative surface potential over the mutation site, not an H-bond network.

**Affinity calibration (Tables 1–2).** The parent JDI TCR binds G12D at **63 µM** with **no detectable binding** to wild-type, and is also clean against G12C, ERAS, RRAS, MRAS, DIRAS2 and RAP1A. As affinity was matured, wild-type binding *appeared*: NB → 25.8 µM → 3.0 µM → 425 nM. **Micromolar affinity is the target, not nanomolar** — over-optimising affinity degrades the objective that is ranked first.

---

# Route A (primary) — conformational selection on the open state

## Rationale

Design a binder that engages an Asp6 projecting out of the groove, using the experimentally observed **open** backbone of the free complex.

The mechanistic claim is conformational selection, not induced fit: the open state need only be *populated* in the free ensemble, not dominant. The penalty on apparent affinity is −RT ln *p* — about 1.4 kcal/mol at 10% population, 2.7 at 1%. A binder that engages the open state shifts the population toward it.

Three things support this:

1. **The disorder in 7OW4 is evidence for the premise.** No density at residues 5–6 in any of four copies means the free G12D peptide centre is genuinely sampling multiple conformations. Conformational selection requires a populated alternative state; disorder is what that looks like in a diffraction experiment.
2. **The discriminator is Asp-versus-nothing.** A buried salt bridge plus H-bonds to a carboxylate is worth 2–5 kcal/mol, and wild-type presents a hydrogen in its place. This is the most favourable substitution pair available, and far more designable than reproducing a 13–16 Å Coulombic coupling.
3. **Negative design comes built in.** A pocket presenting a cation and backbone donors sized for a carboxylate is, against wild-type, a desolvated unsatisfied cation — actively destabilising, not merely neutral. The same feature serves objectives 1 and 3. Objective 2 follows too: apo HLA has an empty groove and nothing for the reader to grab.

**Known headwinds.** The groove rewards burial — only the G12D peptide forms a network with the HLA F-pocket in the bound state, so the closed conformation is *more* stabilised for the mutant than for WT. And the entropy is paid twice: ordering a mobile region, in a non-preferred conformation. The TCR paid TΔS° = −8.13 kcal/mol across ~990 Å² buried; a 60–70 residue miniprotein buries ~800–900 Å², so there is margin but little.

**Expected outcome:** a binder with a real mechanism and a genuine chance at objective 1. Not a high-probability binder.

## Steps

1. **Build the target.** 7OW3 copy 4 — HLA chain **J**, β2m chain **K**, peptide chain **L**. Use L rather than C: φ = −86.4° is more comfortable for a non-Gly residue, and chain L is better ordered throughout (residue 6 mean B 48.9 vs 57.9 Å²). Mutate Gly6 → Asp in the PyMOL mutagenesis wizard; keep the three canonical χ1 wells as separate models.
2. **Protonate and relax.** PDBFixer for missing atoms and hydrogens at pH 7.4; restrained minimisation in OpenMM (amber14-all + implicit GBn2), heavy-atom restraints on HLA and β2m, peptide free. Confirm residue 6 φ/ψ has not left the allowed region and the Asp χ1/χ2 remain canonical.
3. **Test whether the open state is populated — the decisive hour.** Unrestrained MD of the pMHC alone, implicit solvent, 3 × 5 ns (~6,000 atoms; well under an hour each on an M2). Track the residue-6 Cβ→R114 CZ distance and Asp6 side-chain SASA over time. Minimisation alone cannot answer this: it always converges to something strain-free, and strain-free ≠ populated.
4. **Select design frames.** Take 2–3 frames spanning distinct Asp rotamers with side-chain SASA > ~25 Å². Designing against several rotamers hedges the single largest avoidable error — aiming a pocket at empty space.
5. **Generate backbones.** RFdiffusion on Colab, binder mode, target = HLA α1/α2 (chain J res 1–180) + peptide chain L. Hotspots on peptide residues 4–7 plus the Asp. ~200 backbones per frame, `noise_scale_ca=0.5`, `noise_scale_frame=0.5`. Checkpoint to Drive.
6. **Design sequences.** ProteinMPNN, `--sampling_temp 0.1 --omit_AAs C`, 8 sequences per backbone.
7. **Mechanism audit (the primary filter).** For each candidate: confirm the designed reader group is within 3.2 Å of the Asp6 carboxylate and buried (freesasa); then thread the WT peptide into the same model and confirm the reader is left unsatisfied and desolvated. A design that passes has a mechanism; one that fails has a correlation. This filter is cheap, local, and tests the thing being claimed.
8. **Affinity proxy.** PRODIGY on the G12D complex and on the WT-threaded complex; rank by the difference.
9. **Counter-screen.** Co-fold finalists against all three targets with identical settings; rank by ipSAE(G12D) − max(ipSAE(WT), ipSAE(apo)).
10. **Novelty, screens, package.** Foldseek + MMseqs2; the three named screens (pLDDT, ipSAE, mean Kyte–Doolittle hydropathy).

**Route-specific risk — read before trusting step 9.** Co-folding models are trained on the PDB, where this peptide's centre is either absent or closed. They may refuse to reproduce the open state and return poor ipSAE for reasons unrelated to design quality. The selection prompt computes ipSAE itself, with one co-folder at fixed settings, so the screening metric is structurally biased against this mechanism. Mitigation: state this explicitly in the methods; rely on step 7 as primary evidence; and place the best-argued design in row 1, which the prompt protects (*"a failed screen … never removes a participant's ranked top pick"*).

## Tools

PyMOL (open-source) · PDBFixer + OpenMM · RFdiffusion (Colab T4) · ProteinMPNN (local CPU) · freesasa · PRODIGY · Boltz-2 / Chai-1 (Colab) · AlphaFold Server (finalists, full complex) · ipSAE · Foldseek, MMseqs2

---

# Route B (alternative) — closed state, indirect electrostatics

Run only if Route A's step 3 says the open state is not populated.

## Rationale

Design against 7OW6 with the TCR stripped, reproducing Poole's validated mechanism: a paratope whose surface is electrostatically complementary to the more negative potential that the buried Asp6 projects upward.

**What recommends it:** the mechanism is experimentally demonstrated on this exact target, with >4000× discrimination and MMPBSA agreement (12.6 vs 10.8 kcal/mol experimental). The HLA α1/α2 surface is invariant across all five structures, so a binder reading mostly that surface plus the peptide backbone is far more likely to fold and bind.

**What it costs.** Three things to be clear about:

- 7OW6's peptide conformation exists *because a TCR is bound*. This route does not escape the conformational-selection argument; it makes the same bet on a different member of the ensemble.
- The discriminating signal acts from 13–16 Å. That is not deliberately designable with free tools, and not verifiable with them either.
- The invariant HLA surface is identical in all three assay targets, so the default outcome is a binder that hits all of them — failing objectives 1 and 2.

So this route plays for a different prize: **objective-3 honourable mention and getting tested**, not for the win. That is a legitimate bet given the goal hierarchy, but it should not be mistaken for a safer path to the same outcome.

## Steps

1. **Build the target.** 7OW6, keep chains A (HLA), B (β2m), C (peptide with real Asp6); delete chains D and E (TCR). Light restrained minimisation only — the coordinates are experimental.
2. **Generate backbones.** RFdiffusion hotspotted over peptide residues 4–7 and the flanking groove rim (HLA ~65–72 and ~150–158). **Do not hotspot D6 itself** — it is 2.95 Å from R114 and unreachable; aiming at it produces clashing designs.
3. **Design sequences with charge intent.** ProteinMPNN as above, then filter for net-positive formal charge among binder residues within ~12 Å of peptide residue 6. This is the paper's mechanism expressed as a selectable property.
4. **Electrostatic check.** PDB2PQR + APBS on each finalist; confirm a positive potential patch facing the peptide centre. This is the same analysis Poole et al. used (their Fig. 3e,f, ±2 eV) and it is cheap on CPU.
5. **Affinity proxy and counter-screen.** As Route A steps 8–9. Note the caveat below: PRODIGY will be nearly blind to this route's discriminator.
6. **Novelty, screens, package.** As Route A step 10.

## Tools

Same as Route A, minus the MD, plus **PDB2PQR + APBS** for the electrostatic surface.

---

## Decision rule

Set the branch before seeing the result:

- **Open state persists over a few ns** (Cβ→R114 stays > ~9 Å, Asp6 SASA stays > ~20 Å²) → Route A, full allocation.
- **Residue 6 collapses into the groove within ~1 ns** → the population penalty is likely 3–4+ kcal/mol and Route A will bind nothing. Switch to Route B.
- **Mixed / it samples both** → Route A primary, with 2–3 Route B designs in rows 8–10 as a hedge.

Prepare both targets up front (~30 min). The pipelines are otherwise identical, so keep Route B *trailing* rather than parallel — abandonable at any point without having cost the primary.

---

## Affinity prediction — PRODIGY

Xue et al., *Bioinformatics* 2016, [10.1093/bioinformatics/btw514](https://doi.org/10.1093/bioinformatics/btw514). Available locally as `prodigy-prot` (pip) or via the HADDOCK web server.

PRODIGY predicts ΔG from counts of interfacial contacts by residue class plus non-interacting surface area. Two consequences worth recording in the methods:

- **It suits Route A.** Replacing Asp6 with Gly removes side-chain contacts outright, which PRODIGY registers directly.
- **It is near-blind to Route B's discriminator.** If the binder makes no contact with position 6, the contact map barely changes between G12D and WT, so the predicted ΔΔG will be close to zero regardless of whether the design works. For Route B, the APBS electrostatic check is the on-mechanism evidence and PRODIGY only bounds overall affinity.

Trained on natural complexes with ~1.5 kcal/mol error, so use it for **ranking within a congeneric series**, never as an absolute ΔG. Report it as one of several signals; the selection prompt is explicit that no single in silico metric decides.

---

## Day plan (3–8 Nov, ≤2 h/day)

| Day | Attention | Unattended |
|---|---|---|
| **Mon 3** | Fix numpy. Build both targets. Mutate G6→D on 7OW3 chain L, relax. | — |
| **Tue 4** | Launch the MD; read the result; **apply the decision rule**. Set up Colab. | RFdiffusion overnight |
| **Wed 5** | ProteinMPNN; mechanism audit (step 7); PRODIGY triage. | Boltz-2 overnight |
| **Thu 6** | ΔipSAE triage. AlphaFold Server on finalists (full complex, 3 targets each). | AF Server |
| **Fri 7** | Novelty check; compute the three screens; pick and order the 20 rows. | — |
| **Sat 8** | Methods document, package, **submit with hours to spare**. | — |

---

## Submission package

- **CSV**, row order = ranking: `name`, `sequence`, `molecule_class`, plus pLDDT, ipSAE, mean Kyte–Doolittle hydropathy, ΔipSAE, PRODIGY ΔG and ΔΔG, declared epitope, route label.
- **State the epitope per design**, naming contacted residues including the position-6 contact. Unstated epitopes get predicted and read from a model instead.
- **Attach everything as readable text.** Links are never followed — *"what you have no readable copy of counts as not shared."* Export complexes as PDB/mmCIF, not `.pse` (the existing `pdb_models.pse` is unreadable to the selector).
- **Methods document:** the structural analysis above (open vs closed, with the measurements), the conformational-selection argument and the 7OW4 disorder evidence, the MD result and the decision rule as applied, the mechanism audit as primary evidence, the ipSAE bias caveat stated openly, and what could not be done (no GPU, no ProteinTyper output, co-folder limits). Plain description only — anything resembling an instruction to the reader risks being flagged as prompt injection, which sets aside the whole account.

## Open questions

1. **The open-state population is unmeasured.** The MD in Route A step 3 is a short implicit-solvent proxy, not a free-energy calculation. It can falsify the premise cheaply; it cannot confirm it quantitatively. Say so.
2. **The open backbone is wild-type.** 7OW3 carries Gly6; 7OW4's worse disorder suggests the G12D free peptide is, if anything, more mobile there. "The open conformation" is an ensemble average being treated as a state.
3. **Novelty credit is in target preparation, not the algorithm.** RFdiffusion + ProteinMPNN is a crowded method family (25% / 93-design cap per track). What this submission earns is design novelty (level 3/4) plus a checkable mechanism in the write-up.
4. **ipSAE DOI** ([10.1101/2025.02.10.637595](https://doi.org/10.1101/2025.02.10.637595)) is the citation I am least sure of; the GitHub repo is authoritative.
5. **AlphaFold Server limits and terms** change — confirm the per-day quota and non-commercial terms before building a plan on them.
6. **Prior work worth reading:** Sim et al., *PNAS* 2020, [10.1073/pnas.1921964117](https://doi.org/10.1073/pnas.1921964117) — TCRs against this same epitope, and the closest prior art on how receptors achieve mutation specificity here.
