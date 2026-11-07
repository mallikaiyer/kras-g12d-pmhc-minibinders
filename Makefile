# Pipeline driver. Every stage is a script under scripts/, numbered in the
# order it must run. Run `make help` for the list.
#
# Assumes the conda env from environment.yml is active:
#     conda activate kras-pmhc

PYTHON ?= python
SCRIPTS := scripts

.PHONY: help install fetch survey rotamers targets relax test lint clean

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

test:
	$(PYTHON) -m pytest

lint:
	ruff check .

# Deliberately does not touch data/external: re-downloading is slow and the
# manifest exists precisely so the same files can be reused.
clean:
	rm -rf results/survey
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
