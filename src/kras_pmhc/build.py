"""Building side-chain atoms onto a backbone.

Used both to answer "which way could an Asp at position 6 actually point?" and,
in the target-preparation stage, to place the mutant side chain itself.

Atoms are placed by internal coordinates (bond length, bond angle, torsion)
using the natural-extension-reference-frame construction, which is the same
arithmetic PyMOL's mutagenesis wizard and every rotamer library use.
"""

from __future__ import annotations

import numpy as np
from biotite.structure import AtomArray

from .geometry import virtual_cb

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


#: Atoms added when turning a glycine into an aspartate, in the order a PDB
#: file conventionally lists them.
ASP_ADDED_ATOMS = ("CB", "CG", "OD1", "OD2")


def mutate_gly_to_asp(
    array: AtomArray,
    chain_id: str,
    res_id: int,
    chi1: float,
    chi2: float,
) -> AtomArray:
    """Replace a glycine with an aspartate at a chosen rotamer.

    The backbone is passed through untouched: Cbeta is reconstructed from the
    existing N-CA-C frame, so the mutation adds four atoms and relabels the
    residue without moving anything that was observed experimentally. That
    matters because the whole argument for this target is that the backbone is
    an experimental observation rather than a model.

    Args:
        array: structure containing the glycine. Not modified.
        chain_id, res_id: which residue to mutate.
        chi1: N-CA-CB-CG torsion in degrees.
        chi2: CA-CB-CG-OD1 torsion in degrees.

    Returns:
        A new array with CB, CG, OD1 and OD2 inserted directly after the
        residue's existing atoms, so that output files stay grouped by residue.

    Raises:
        ValueError: if the residue is absent, is not a glycine, or lacks a
            complete backbone.
    """
    mask = (array.chain_id == chain_id) & (array.res_id == res_id)
    target = array[mask]
    if len(target) == 0:
        raise ValueError(f"chain {chain_id} residue {res_id} is not modelled")
    observed = str(target.res_name[0])
    if observed != "GLY":
        raise ValueError(
            f"chain {chain_id} residue {res_id} is {observed}, expected GLY"
        )

    coords = {}
    for name in ("N", "CA", "C"):
        hit = target.coord[target.atom_name == name]
        if len(hit) == 0:
            raise ValueError(f"chain {chain_id} residue {res_id} has no {name} atom")
        coords[name] = np.asarray(hit[0], dtype=float)

    cb = virtual_cb(coords["N"], coords["CA"], coords["C"])
    built = {"CB": cb, **build_asp_sidechain(coords["N"], coords["CA"], cb, chi1, chi2)}

    work = array.copy()
    work.res_name[mask] = "ASP"

    # Duplicating an existing atom guarantees the new atoms carry exactly the
    # same annotation categories as the rest of the array, including optional
    # ones such as b_factor and occupancy. Building a bare AtomArray instead
    # would drop them and make concatenation fail.
    indices = np.where(mask)[0]
    extra = work[np.repeat(indices[0], len(ASP_ADDED_ATOMS))]
    for i, name in enumerate(ASP_ADDED_ATOMS):
        extra.coord[i] = built[name]
        extra.atom_name[i] = name
        extra.element[i] = name[0]
    if "b_factor" in extra.get_annotation_categories():
        # Built atoms inherit the local disorder of the residue they extend,
        # which is more honest than claiming a B-factor of zero.
        extra.b_factor[:] = float(np.mean(target.b_factor))

    insert_at = int(indices[-1]) + 1
    return work[:insert_at] + extra + work[insert_at:]
