# Colab runbook — stage 05 backbone generation

RFdiffusion needs a CUDA GPU; this project's machine is a MacBook Air M2. Stage 05
is therefore the one stage that is not fully scripted on this machine, and this
document records exactly what is done by hand so the run is reproducible and
auditable.

Everything either side of the notebook *is* scripted: `make stage` puts the target
where Colab can read it, and `make backbones` reads the results back and scores them.
The manual part is confined to opening a browser, selecting a GPU, and running cells.

## Why Google Drive rather than an upload cell

Google Drive for Desktop is mounted on this machine, so staging the target is a plain
file copy and needs no OAuth client, no API credentials and nothing embedded in the
notebook. The same folder carries results back: the notebook writes each backbone to
Drive, Drive syncs it here, and the ingest script reads it directly. No zip, no
manual download, and the notebook stays a small text file.

`scripts/05a_stage_to_drive.py` writes `TARGETS.tsv` alongside the target, recording
its SHA-256. Notebook cell 3 re-checks that checksum. This is what ties generated
backbones to a specific stage-04 output: a wrong or half-synced file fails loudly
instead of producing plausible backbones against the wrong geometry.

## Procedure

### 1. Stage the target (scripted, on the Mac)

```bash
conda activate kras-pmhc
make stage
```

Confirm it printed `copied:` (or `unchanged:`) and `wrote: TARGETS.tsv`.

### 2. Wait for Drive to sync — manual

Watch the Google Drive icon in the macOS menu bar until it reports sync complete.
The target is 124 KB, so this is seconds, but cell 3 will fail on a partial file.

### 3. Open the notebook and attach a GPU — manual

Open `notebooks/05_rfdiffusion_binder.ipynb` in Colab
(File → Upload notebook, or via GitHub once the repo is published).

**Runtime → Change runtime type → T4 GPU → Save.**

### 4. Run cells 1–5 — manual

| Cell | What it does | Expect |
|---|---|---|
| 1 | Parameters | prints contigs, hotspots, noise, cutoff |
| 2 | Install RFdiffusion + ColabDesign, download checkpoints | ~3 min; ends with `nvidia-smi` naming the GPU |
| 3 | Mount Drive, load target, verify SHA-256 | `sha256 ... verified against TARGETS.tsv` |
| 4 | Define the generation driver | prints the contig string and `resuming at batch 0` |
| 5 | **Pilot: 8 designs** | the full `run_inference.py` command, then 8 Cα traces |

Record below: the GPU allocated, the session start time, and the per-design time
printed by cell 5.

If cell 5 fails with a hydra error naming `inference.seed`, set `USE_SEED = False` in
cell 1 and re-run cells 1, 4, 5. Note that outputs are then unseeded.

### 5. Ingest and read the pilot gate (scripted, on the Mac)

```bash
make backbones
```

The gate is applied here rather than in the notebook because it uses the Levy (2010)
rSASA classification from `kras_pmhc.interface`, and having two implementations of
the project's interface definition would be worse than a round trip.

A backbone passes when **Asp6 is an interface residue** (support, rim or core) **and
the peptide carries ≥20% of the target-side interface area**. Proceed only if **≥4
of 8** pass.

The script exits non-zero on a hard failure — a moved target, a wrong binder length,
a steric overlap, or a one-sided interface. Those are bugs or wrong inputs, not design
failures; do not continue past them.

**If fewer than 4 pass**, the binder is landing on the α1/α2 helix tops rather than
the peptide. That is the predicted failure mode: the helix tops are a far larger
hydrophobic surface than a ten-residue peptide in a groove. Remedies, one pilot each:

1. widen `HOTSPOTS` to `"C4,C5,C6,C7,C8"`;
2. then set `use_beta_model` — switch `CHECKPOINT` to `Complex_beta_ckpt.pt` — for
   less helical topologies.

### 6. Production — manual

Only if the gate passed. Set `CUTOFF_UTC` in cell 1 (Colab runs UTC; the competition
deadline is Mon 9 Nov 06:59 EST = 11:59 UTC), re-run cell 1, then run cell 6. It
generates in batches of 8 and checkpoints each batch to Drive, so a preemption loses
at most one batch and cell 6 can simply be re-run to resume.

### 7. Re-ingest — scripted

```bash
make backbones
```

## Session log

Fill in as the run proceeds. This is the record of what actually happened, not what
was planned.

| Field | Value |
|---|---|
| Date |  |
| GPU allocated |  |
| Session start (UTC) |  |
| Pilot: per-design time |  |
| Pilot: passing / 8 |  |
| Hotspots used |  |
| Checkpoint used |  |
| `inference.seed` accepted? |  |
| `CUTOFF_UTC` set to |  |
| Preemptions (count, times) |  |
| Total backbones generated |  |
| Notes |  |

## Known deviations from the plan

Record any here, with the reason. An empty section means none.
