# Pipeline driver. Every stage is a script under scripts/, numbered in the
# order it must run. Run `make help` for the list.
#
# Assumes the conda env from environment.yml is active:
#     conda activate kras-pmhc

PYTHON ?= python
SCRIPTS := scripts

.PHONY: help install fetch survey rotamers targets relax stage backbones nim nim-key nim-smoke rotamer-sensitivity sequences select stage-cofolding test lint clean

# Where scripts/05a_stage_to_drive.py puts the target for Colab to read, and
# where scripts/05_ingest_backbones.py reads the results back from. Override on
# the command line if Drive is mounted elsewhere:
#     make stage DRIVE_DIR="/path/to/My Drive/kras_pmhc_stage05"
DRIVE_DIR ?= $(HOME)/Library/CloudStorage/GoogleDrive-you@example.com/My Drive/kras_pmhc_stage05
DRIVE_DIR_07 ?= $(HOME)/Library/CloudStorage/GoogleDrive-you@example.com/My Drive/kras_pmhc_stage07

help:
	@echo "Setup"
	@echo "  make install   install this package into the active env (editable)"
	@echo "  make fetch     download the five PDB entries and write the manifest"
	@echo ""
	@echo "Analysis"
	@echo "  make survey    structural survey -> results/survey/"
	@echo "  make rotamers  Asp rotamer scan -> results/rotamer_scan/"
	@echo "  make targets   build design targets -> data/targets/"
	@echo "  make relax     protonate + minimise targets -> data/targets_relaxed/"
	@echo ""
	@echo "Generation (stage 05; needs a GPU, so the middle step runs in Colab)"
	@echo "  make stage     copy the target to Drive for notebooks/05_rfdiffusion_binder.ipynb"
	@echo "  make backbones ingest generated backbones -> data/backbones/"
	@echo "                 (add RELABEL=1 for raw output from the Kaggle route)"
	@echo "  make rotamer-sensitivity  re-score designs against an alternative Asp6 rotamer"
	@echo "  make nim-key   store the NVIDIA API key (hidden prompt)"
	@echo "  make nim-smoke one design via the NVIDIA NIM API (proves the request)"
	@echo "  make nim       10 designs via the NVIDIA NIM API (one call each)"
	@echo ""
	@echo ""
	@echo "Sequence design (stage 06)"
	@echo "  make sequences ProteinMPNN sequences for the designs that passed -> data/sequences/"
	@echo "  make select    drop low-complexity and lone-cysteine sequences -> SELECTED.tsv"
	@echo ""
	@echo "Cofolding (stage 07; needs a GPU, so the middle step runs in Colab)"
	@echo "  make stage-cofolding  write Boltz-2 inputs to Drive for notebooks/07_*.ipynb"
	@echo ""
	@echo "Quality"
	@echo "  make test      run the test suite"
	@echo "  make lint      run ruff"
	@echo "  make clean     remove generated results (keeps fetched data)"

install:
	$(PYTHON) -m pip install -e ".[dev]"

fetch:
	$(PYTHON) $(SCRIPTS)/00_fetch_structures.py

survey:
	$(PYTHON) $(SCRIPTS)/01_structure_survey.py

rotamers:
	$(PYTHON) $(SCRIPTS)/02_asp_rotamer_scan.py

targets:
	$(PYTHON) $(SCRIPTS)/03_build_targets.py

relax:
	$(PYTHON) $(SCRIPTS)/04_relax_targets.py

stage:
	$(PYTHON) $(SCRIPTS)/05a_stage_to_drive.py --drive-dir "$(DRIVE_DIR)"

# Raw RFdiffusion output is one continuously renumbered chain. The Kaggle
# route writes it verbatim, so it needs relabelling at ingest; the Colab route
# relabelled in the notebook and does not. Pass RELABEL=1 for the former:
#     make backbones DRIVE_DIR=<unzipped>/stage05 RELABEL=1
RELABEL ?=
CONTIGS ?= A1-180/0 C1-10/0 60-75

backbones:
	$(PYTHON) $(SCRIPTS)/05_ingest_backbones.py --drive-dir "$(DRIVE_DIR)" \
	  $(if $(RELABEL),--relabel-contigs "$(CONTIGS)",)

# Fallback transport for stage 05 when no GPU is available. Needs an API key in
# .nim_api_key (gitignored). NIM_DIR is where backbones and run_log.tsv land, in
# the same layout the Colab notebook wrote, so `backbones` can ingest it:
#     make backbones DRIVE_DIR=data/nim_stage05
NIM_DIR ?= data/nim_stage05
NIM_N ?= 10

# Quantifies how far the stage-05 gate depends on the constructed Asp6 rotamer.
# Writes results/rotamer_sensitivity.tsv; changes no decision.
# Best N sequences per design, three predictions each (g12d / wt / apo).
PER_DESIGN ?= 2

stage-cofolding:
	$(PYTHON) $(SCRIPTS)/07a_stage_cofolding.py --drive-dir "$(DRIVE_DIR_07)" \
	  --per-design $(PER_DESIGN) --force

select:
	$(PYTHON) $(SCRIPTS)/06b_select_sequences.py

sequences:
	$(PYTHON) $(SCRIPTS)/06_mpnn_sequences.py

rotamer-sensitivity:
	$(PYTHON) $(SCRIPTS)/05c_rotamer_sensitivity.py

nim-key:
	$(PYTHON) $(SCRIPTS)/nim_set_key.py

nim-smoke:
	$(PYTHON) $(SCRIPTS)/05b_nim_generate.py --n 1 --out-dir "$(NIM_DIR)"

nim:
	$(PYTHON) $(SCRIPTS)/05b_nim_generate.py --n $(NIM_N) --out-dir "$(NIM_DIR)"

test:
	$(PYTHON) -m pytest

lint:
	ruff check .

# Deliberately does not touch data/external: re-downloading is slow and the
# manifest exists precisely so the same files can be reused.
clean:
	rm -rf results/survey
	find . -name __pycache__ -type d -prune -exec rm -rf {} +

