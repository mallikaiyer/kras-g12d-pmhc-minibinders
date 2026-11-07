"""Building side-chain atoms onto a backbone.

Used both to answer "which way could an Asp at position 6 actually point?" and,
in the target-preparation stage, to place the mutant side chain itself.

Atoms are placed by internal coordinates (bond length, bond angle, torsion)
using the natural-extension-reference-frame construction, which is the same
arithmetic PyMOL's mutagenesis wizard and every rotamer library use.
"""

from __future__ import annotations

import numpy as np

# Ideal internal coordinates for an aspartate side chain, from standard
# residue topology (CHARMM/AMBER agree to within a few hundredths of an
# angstrom and a degree, which is far below the uncertainty that matters here).
CA_CB = 1.530
CB_CG = 1.522
CG_OD = 1.250
ANGLE_CA_CB_CG = 114.0
ANGLE_CB_CG_OD = 118.4

#: Canonical Asp chi1 rotamer wells (Dunbrack backbone-independent library).
#: The first is by far the most populated in folded proteins.
ASP_CHI1_WELLS = (-70.0, -177.0, 62.0)

#: The carboxylate is planar and its two oxygens are chemically equivalent, so
#: chi2 repeats every 180 degrees. Scanning beyond that duplicates geometry.
ASP_CHI2_PERIOD = 180.0


def place_atom(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    bond: float,
    angle_deg: float,
    torsion_deg: float,
) -> np.ndarray:
    """Place a fourth atom D given three existing atoms A-B-C.

    Args:
        a, b, c: coordinates of the three reference atoms, shape ``(3,)``.
        bond: C-D bond length in angstroms.
        angle_deg: B-C-D bond angle in degrees.
        torsion_deg: A-B-C-D torsion in degrees.

    Returns:
        Coordinates of D, shape ``(3,)``.
    """
    a, b, c = (np.asarray(v, dtype=float) for v in (a, b, c))
    angle = np.deg2rad(angle_deg)
    torsion = np.deg2rad(torsion_deg)

    # Orthonormal frame at C: bc along the previous bond, n normal to the
    # A-B-C plane, and their cross product completing a right-handed basis.
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    norm_n = np.linalg.norm(n)
    if norm_n < 1e-8:
        raise ValueError("reference atoms A, B, C are collinear; torsion is undefined")
    n /= norm_n
    m = np.cross(n, bc)

    # D in the local frame, then rotated into place.
    d_local = np.array(
        [
            -bond * np.cos(angle),
            bond * np.sin(angle) * np.cos(torsion),
            bond * np.sin(angle) * np.sin(torsion),
        ]
    )
    return c + d_local[0] * bc + d_local[1] * m + d_local[2] * n


def build_asp_sidechain(
    n: np.ndarray,
    ca: np.ndarray,
    cb: np.ndarray,
    chi1: float,
    chi2: float,
) -> dict[str, np.ndarray]:
    """Construct CG, OD1 and OD2 for an aspartate at a given rotamer.

    Args:
        n, ca, cb: backbone nitrogen, alpha carbon, and beta carbon. ``cb`` may
            be a reconstructed virtual Cbeta, which is how a glycine backbone is
            used as a template for a mutant side chain.
        chi1: N-CA-CB-CG torsion, degrees.
        chi2: CA-CB-CG-OD1 torsion, degrees.

    Returns:
        Mapping of atom name to coordinates.
    """
    cg = place_atom(n, ca, cb, CB_CG, ANGLE_CA_CB_CG, chi1)
    od1 = place_atom(ca, cb, cg, CG_OD, ANGLE_CB_CG_OD, chi2)
    # The second carboxylate oxygen is the first reflected through the plane,
    # i.e. 180 degrees further round the CB-CG bond.
    od2 = place_atom(ca, cb, cg, CG_OD, ANGLE_CB_CG_OD, chi2 + 180.0)
    return {"CG": cg, "OD1": od1, "OD2": od2}
