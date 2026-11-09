"""Checks on the stage-05 notebooks.

A notebook cannot be imported or exercised locally - it needs a GPU and either
a Drive mount or a Kaggle dataset - so a syntax error in it would only surface
after the setup cell has spent minutes installing RFdiffusion, on a free GPU
session that may not be granted again. These tests are cheap insurance against
wasting one, and the morning of 9 Nov 2026 is what they are insurance against.

Ruff does not lint ``notebooks/`` (see pyproject.toml): the Colab setup cell is
reproduced verbatim from upstream and must stay that way.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

NOTEBOOKS = Path(__file__).resolve().parents[1] / "notebooks"
COLAB = NOTEBOOKS / "05_rfdiffusion_binder.ipynb"
KAGGLE = NOTEBOOKS / "05_kaggle_rfdiffusion_binder.ipynb"

#: Line prefixes that are IPython syntax rather than Python, and so must be
#: stripped before the cell can be parsed.
IPYTHON_PREFIXES = ("%", "!")


def _load(path: Path) -> dict:
    if not path.exists():
        pytest.fail(f"{path.name} is missing")
    return json.loads(path.read_text())


def _source(cell) -> str:
    source = cell["source"]
    return source if isinstance(source, str) else "".join(source)


def _code_cells(notebook):
    return [c for c in notebook["cells"] if c["cell_type"] == "code"]


def _all_code(notebook) -> str:
    return "\n".join(_source(c) for c in _code_cells(notebook))


@pytest.fixture(scope="module")
def colab():
    return _load(COLAB)


@pytest.fixture(scope="module")
def kaggle():
    return _load(KAGGLE)


# --- both notebooks -------------------------------------------------------


@pytest.mark.parametrize("path", [COLAB, KAGGLE], ids=["colab", "kaggle"])
def test_notebook_is_valid_nbformat_and_asks_for_a_gpu(path):
    notebook = _load(path)
    assert notebook["nbformat"] == 4
    assert notebook["cells"], "notebook has no cells"
    assert notebook["metadata"]["accelerator"] == "GPU", (
        "the notebook must request a GPU runtime"
    )


@pytest.mark.parametrize("path", [COLAB, KAGGLE], ids=["colab", "kaggle"])
def test_every_code_cell_parses(path):
    """The check that protects a GPU session from a typo."""
    notebook = _load(path)
    for index, cell in enumerate(_code_cells(notebook)):
        lines = [
            line
            for line in _source(cell).split("\n")
            if not line.strip().startswith(IPYTHON_PREFIXES)
        ]
        try:
            ast.parse("\n".join(lines))
        except SyntaxError as error:
            pytest.fail(f"{path.name} code cell {index} does not parse: {error}")


@pytest.mark.parametrize("path", [COLAB, KAGGLE], ids=["colab", "kaggle"])
def test_the_parameters_cell_pins_the_agreed_values(path):
    """Guards the decisions in docs/design_plan.md against silent drift.

    These are not arbitrary: 0.5 is the RFdiffusion README's recommendation for
    protein-protein interfaces, the hotspots are centred on the mutation, and
    the contigs are the cropped stage-04 target. A change to any of them is a
    change to the experiment and should have to be made deliberately. Both
    routes must agree on them, or the two sets of backbones are not comparable.
    """
    source = _source(_code_cells(_load(path))[0])
    for expected in (
        'HOTSPOTS = "C5,C6,C7"',
        "NOISE_SCALE = 0.5",
        'BINDER_LENGTH',
        'CHECKPOINT',
        '"A1-180/0", "C1-10/0"',
        "7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb",
        'CUTOFF_UTC = "2026-11-09T20:00:00"',
    ):
        assert expected in source, f"{path.name}: parameters cell no longer sets {expected}"


@pytest.mark.parametrize("path", [COLAB, KAGGLE], ids=["colab", "kaggle"])
def test_the_target_checksum_is_enforced(path):
    """A wrong or half-synced target must fail loudly, not silently."""
    source = _all_code(_load(path))
    assert "sha256" in source.lower()
    assert "mismatch" in source


@pytest.mark.parametrize("path", [COLAB, KAGGLE], ids=["colab", "kaggle"])
def test_the_generation_driver_passes_every_flag_explicitly(path):
    """The printed command must be a complete record of what ran.

    If generation went through a helper that assembles the flags, the noise
    scale could not be set and the printed command would no longer account for
    the coordinates.
    """
    source = _all_code(_load(path))
    for flag in (
        "denoiser.noise_scale_ca",
        "denoiser.noise_scale_frame",
        "contigmap.contigs",
        "ppi.hotspot_res",
        "inference.ckpt_override_path",
        "diffuser.T",
    ):
        assert flag in source, f"{path.name}: {flag} is not set anywhere"
    assert "run_diffusion(" not in source, (
        "run_diffusion hides the hydra flags; the driver must call "
        "run_inference.py directly"
    )


# --- Kaggle-specific ------------------------------------------------------


def test_the_kaggle_notebook_refuses_python_313(kaggle):
    """The failure that cost 9 Nov 2026: on 3.13 pip resolves dgl to 0.1.3.

    The notebook must build a 3.12 environment rather than trust the image.
    """
    source = _all_code(kaggle)
    assert "sys.version_info[:2] <= (3, 12)" in source
    assert "uv venv --python 3.12" in source


def test_the_kaggle_notebook_proves_dgl_reaches_the_gpu_before_generating(kaggle):
    """Importing dgl is not enough: a CPU-only build fails at ``g.to('cuda')``.

    The probe must run before the pilot, so a wrong environment costs seconds
    instead of a batch.
    """
    source = _all_code(kaggle)
    assert "g.to('cuda')" in source
    assert "do not run the pilot" in source


def test_the_kaggle_notebook_does_not_relabel_chains_itself(kaggle):
    """One tested implementation of this, not one per transport.

    Raw RFdiffusion output is relabelled on the Mac by
    ``kras_pmhc.nim.relabel_chains`` via ``--relabel-contigs``.
    """
    source = _all_code(kaggle)
    assert "def fix_chains" not in source
    assert "relabel" in source.lower()


def test_the_kaggle_archive_carries_the_run_log(kaggle):
    """The ingest script refuses to run without run_log.tsv."""
    source = _all_code(kaggle)
    assert "run_log.tsv" in source
    assert "make_archive" in source
