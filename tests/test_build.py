"""Tests for side-chain construction.

The side-chain builder is validated the same way as the virtual Cbeta: by
rebuilding a real aspartate from its own torsion angles and comparing with the
deposited coordinates. If this drifts, the rotamer scan that justifies Route A
is measuring the wrong geometry.
"""

from __future__ import annotations

import biotite.structure as struc
import numpy as np
import pytest

from kras_pmhc import structures
from kras_pmhc.build import build_asp_sidechain, place_atom
from kras_pmhc.constants import MUTANT_POSITION


def test_rebuilds_an_observed_aspartate(ow6):
    """Rebuild Asp6 of 7OW6 from its measured chi1/chi2."""
    get = lambda name: structures.atom_coord(ow6, "C", MUTANT_POSITION, name)  # noqa: E731
    n, ca, cb = get("N"), get("CA"), get("CB")

    chi1 = float(np.rad2deg(struc.dihedral(n, ca, cb, get("CG"))))
    chi2 = float(np.rad2deg(struc.dihedral(ca, cb, get("CG"), get("OD1"))))
    # Independently measured in PyMOL as -68.3 / -20.0.
    assert chi1 == pytest.approx(-68.3, abs=0.5)
    assert chi2 == pytest.approx(-20.0, abs=0.5)

    built = build_asp_sidechain(n, ca, cb, chi1, chi2)
    for name in ("CG", "OD1", "OD2"):
        error = float(np.linalg.norm(built[name] - get(name)))
        # Ideal bond lengths and angles against a refined 2.6 A model; the
        # carboxylate oxygens carry the most deviation because the real group
        # is not perfectly planar-symmetric.
        assert error < 0.2, f"{name} rebuilt {error:.3f} A from its observed position"


def test_place_atom_reproduces_requested_internal_coordinates():
    """A placed atom should have the bond, angle and torsion it was asked for."""
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([1.5, 0.0, 0.0])
    c = np.array([2.0, 1.4, 0.0])
    d = place_atom(a, b, c, bond=1.52, angle_deg=114.0, torsion_deg=-70.0)

    assert float(np.linalg.norm(d - c)) == pytest.approx(1.52, abs=1e-6)
    angle = np.rad2deg(struc.angle(b, c, d))
    assert float(angle) == pytest.approx(114.0, abs=1e-4)
    torsion = np.rad2deg(struc.dihedral(a, b, c, d))
    assert float(torsion) == pytest.approx(-70.0, abs=1e-4)


def test_place_atom_rejects_collinear_reference():
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([1.0, 0.0, 0.0])
    c = np.array([2.0, 0.0, 0.0])
    with pytest.raises(ValueError, match="collinear"):
        place_atom(a, b, c, 1.5, 114.0, 60.0)
