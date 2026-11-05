"""Shared fixtures. Structures are loaded once per session - parsing is the
slowest part of the suite and none of the tests mutate the arrays.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from kras_pmhc import structures

PDB_DIR = Path(__file__).resolve().parents[1] / "data" / "external" / "pdb"


def _load_or_skip(pdb_id: str):
    path = PDB_DIR / f"{pdb_id.lower()}.cif"
    if not path.exists():
        pytest.skip(f"{pdb_id} not fetched; run `make fetch`")
    return structures.load(pdb_id)


@pytest.fixture(scope="session")
def ow3():
    """Wild-type peptide, free complex. Four copies, two complete at position 6."""
    return _load_or_skip("7OW3")


@pytest.fixture(scope="session")
def ow4():
    """G12D peptide, free complex. Four copies, none modelled at position 6."""
    return _load_or_skip("7OW4")


@pytest.fixture(scope="session")
def ow5():
    """Wild-type peptide, JDIa41b1 TCR bound."""
    return _load_or_skip("7OW5")


@pytest.fixture(scope="session")
def ow6():
    """G12D peptide, JDIa41b1 TCR bound. The only entry with an ordered Asp6
    whose per-residue energetics are published."""
    return _load_or_skip("7OW6")


@pytest.fixture(scope="session")
def pb2():
    """G12D peptide, parent JDI TCR bound. Two copies."""
    return _load_or_skip("7PB2")
