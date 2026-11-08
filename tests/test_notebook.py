"""Checks on the stage-05 Colab notebook.

A notebook cannot be imported or exercised locally - it needs a GPU and a Drive
mount - so a syntax error in it would only surface after the setup cell has
spent three minutes installing RFdiffusion, on a free GPU session that may not
be granted again. These tests are cheap insurance against wasting one.

Ruff does not lint ``notebooks/`` (see pyproject.toml): the setup cell is
reproduced verbatim from upstream and must stay that way.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest

NOTEBOOK = (
    Path(__file__).resolve().parents[1]
    / "notebooks" / "05_rfdiffusion_binder.ipynb"
)

#: Line prefixes that are IPython syntax rather than Python, and so must be
#: stripped before the cell can be parsed.
IPYTHON_PREFIXES = ("%", "!")


@pytest.fixture(scope="module")
def notebook():
    if not NOTEBOOK.exists():
        pytest.fail(f"{NOTEBOOK.name} is missing")
    return json.loads(NOTEBOOK.read_text())


def _source(cell) -> str:
    source = cell["source"]
    return source if isinstance(source, str) else "".join(source)


def _code_cells(notebook):
    return [c for c in notebook["cells"] if c["cell_type"] == "code"]


def test_notebook_is_valid_nbformat(notebook):
    assert notebook["nbformat"] == 4
    assert notebook["cells"], "notebook has no cells"
    assert notebook["metadata"]["accelerator"] == "GPU", (
        "the notebook must request a GPU runtime"
    )


def test_every_code_cell_parses(notebook):
    """The check that protects a Colab session from a typo."""
    for index, cell in enumerate(_code_cells(notebook)):
        lines = [
            line
            for line in _source(cell).split("\n")
            if not line.strip().startswith(IPYTHON_PREFIXES)
        ]
        try:
            ast.parse("\n".join(lines))
        except SyntaxError as error:
            pytest.fail(f"code cell {index} does not parse: {error}")


def test_the_parameters_cell_pins_the_agreed_values(notebook):
    """Guards the decisions in docs/design_plan.md against silent drift.

    These are not arbitrary: 0.5 is the RFdiffusion README's recommendation for
    protein-protein interfaces, the hotspots are centred on the mutation, and
    the contigs are the cropped stage-04 target. A change to any of them is a
    change to the experiment and should have to be made deliberately.
    """
    source = _source(_code_cells(notebook)[0])
    for expected in (
        'HOTSPOTS = "C5,C6,C7"',
        "NOISE_SCALE = 0.5",
        "DIFFUSER_T  = 50",
        'BINDER_LENGTH  = "60-75"',
        'CHECKPOINT  = "Complex_base_ckpt.pt"',
        '"A1-180/0", "C1-10/0"',
        '7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb',
    ):
        assert expected in source, f"parameters cell no longer sets {expected}"


def test_the_target_checksum_is_enforced(notebook):
    """A wrong or half-synced target must fail loudly, not silently."""
    sources = [_source(c) for c in _code_cells(notebook)]
    assert any("TARGETS.tsv" in s and "sha256" in s for s in sources)
    assert any("checksum mismatch" in s for s in sources)


def test_the_generation_driver_passes_every_flag_explicitly(notebook):
    """run_diffusion is deliberately not used, so the flags must be visible.

    If generation went through ColabDesign's ``run_diffusion`` instead, the
    noise scale could not be set and the printed command would no longer be a
    complete record of what produced the coordinates.
    """
    source = "\n".join(_source(c) for c in _code_cells(notebook))
    for flag in (
        "denoiser.noise_scale_ca",
        "denoiser.noise_scale_frame",
        "contigmap.contigs",
        "ppi.hotspot_res",
        "inference.ckpt_override_path",
        "diffuser.T",
    ):
        assert flag in source, f"{flag} is not set anywhere in the notebook"
    assert "run_diffusion(" not in source, (
        "run_diffusion hides the hydra flags; the driver must call "
        "run_inference.py directly"
    )
