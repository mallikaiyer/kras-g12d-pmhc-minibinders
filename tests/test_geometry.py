"""Tests for the geometry helpers.

The virtual-Cbeta construction is the load-bearing piece: the open-versus-closed
comparison depends on reconstructing a notional side-chain direction for a
glycine and treating it as comparable to an observed aspartate. If that
reconstruction were wrong, the central structural claim of Route A would be
wrong with it, so it is validated against a real Cbeta from a deposited model.
"""

from __future__ import annotations

import numpy as np
import pytest

from kras_pmhc import structures
from kras_pmhc.constants import MUTANT_POSITION
from kras_pmhc.geometry import min_distance, virtual_cb


def test_virtual_cb_reproduces_an_observed_cb(ow6):
    """Reconstruct Asp6's Cbeta in 7OW6 and compare with the deposited atom."""
    n = structures.atom_coord(ow6, "C", MUTANT_POSITION, "N")
    ca = structures.atom_coord(ow6, "C", MUTANT_POSITION, "CA")
    c = structures.atom_coord(ow6, "C", MUTANT_POSITION, "C")
    observed = structures.atom_coord(ow6, "C", MUTANT_POSITION, "CB")
    assert observed is not None, "7OW6 Asp6 should have an observed CB"

    error = float(np.linalg.norm(virtual_cb(n, ca, c) - observed))
    # Ideal-geometry reconstruction against a real 2.6 A model. Agreement is
    # typically ~0.1 A; 0.25 A is a generous ceiling that would still catch a
    # sign error or a transposed coefficient.
    assert error < 0.25, f"virtual CB off by {error:.3f} A"


def test_virtual_cb_rejects_wrong_shape():
    good = np.zeros(3)
    with pytest.raises(ValueError, match="shape"):
        virtual_cb(np.zeros((2, 3)), good, good)


def test_min_distance_simple_case():
    a = np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    b = np.array([[3.0, 4.0, 0.0]])
    assert min_distance(a, b) == pytest.approx(5.0)


def test_min_distance_rejects_empty():
    with pytest.raises(ValueError, match="non-empty"):
        min_distance(np.zeros((0, 3)), np.zeros((1, 3)))
