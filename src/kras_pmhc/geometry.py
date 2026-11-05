"""Geometry helpers not provided by biotite.

Distances, dihedrals and SASA all come from ``biotite.structure``; only the
virtual-Cbeta construction is implemented here, because biotite has no
equivalent and the comparison it enables is central to this project.
"""

from __future__ import annotations

import numpy as np

# Ideal-geometry coefficients for placing Cbeta from the backbone frame, as
# used by trRosetta and AlphaFold. They let glycine - which has no Cbeta - be
# assigned a notional side-chain direction, which is the only way to compare
# the orientation of peptide position 6 between a wild-type Gly structure and
# a mutant Asp structure on equal terms.
#
# Accuracy is checked in tests/test_geometry.py by reconstructing a real Cbeta
# from a deposited structure, which agrees to within ~0.1 A.
_CB_COEFF_NORMAL = -0.58273431
_CB_COEFF_CA_N = 0.56802827
_CB_COEFF_C_CA = -0.54067466


def virtual_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Reconstruct the Cbeta position implied by a backbone N-CA-C frame.

    Args:
        n, ca, c: Cartesian coordinates, shape ``(3,)``, of the backbone
            nitrogen, alpha carbon and carbonyl carbon of one residue.

    Returns:
        Coordinates of the implied Cbeta, shape ``(3,)``.
    """
    n, ca, c = (np.asarray(v, dtype=float) for v in (n, ca, c))
    for name, v in (("n", n), ("ca", ca), ("c", c)):
        if v.shape != (3,):
            raise ValueError(f"{name} must have shape (3,), got {v.shape}")

    b = ca - n
    c_vec = c - ca
    normal = np.cross(b, c_vec)
    return (
        _CB_COEFF_NORMAL * normal
        + _CB_COEFF_CA_N * b
        + _CB_COEFF_C_CA * c_vec
        + ca
    )


def min_distance(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Closest approach between two sets of coordinates, in angstroms.

    Args:
        coords_a, coords_b: arrays of shape ``(n, 3)`` and ``(m, 3)``.

    Raises:
        ValueError: if either set is empty. Returning infinity instead would
            silently poison a contact table with a plausible-looking number.
    """
    a = np.atleast_2d(np.asarray(coords_a, dtype=float))
    b = np.atleast_2d(np.asarray(coords_b, dtype=float))
    if a.size == 0 or b.size == 0:
        raise ValueError("min_distance requires two non-empty coordinate sets")
    # Pairwise distances via broadcasting; these groups are at most a few dozen
    # atoms, so the full matrix is cheap.
    return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).min())
