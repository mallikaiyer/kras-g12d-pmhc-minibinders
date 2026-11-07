"""The four structural analyses the design plan rests on.

Each function returns a list of plain dicts, one per row, so that callers can
write TSV, build a DataFrame or assert in tests without this module knowing
about output formats.

The analyses, and what each decides:

1. ``peptide_density``  - which copies of the peptide have modelled coordinates
   at the central positions. Establishes that 7OW4 has none in any copy, and
   that 7OW3 has two complete copies out of four.
2. ``peptide_dihedrals`` - whether the backbone an Asp would inherit is
   compatible with a non-glycine residue. A positive phi would be
   glycine-only territory and would sink Route A.
3. ``position6_orientation`` - whether position 6 points into the groove or out
   of it, per copy. This is the open-versus-closed measurement.
4. ``mutant_contacts`` - how far the Asp6 carboxylate is from the TCR residues
   that carry the discrimination, and from its HLA partners. Establishes that
   the mechanism is long-range rather than contact-mediated.
"""

from __future__ import annotations

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray

from . import structures
from .constants import (
    ASP_CARBOXYLATE_ATOMS,
    GROOVE_REFERENCE_ATOM,
    HLA_SURVEY_RESIDUES,
    MUTANT_POSITION,
    STRUCTURES,
    TCR_DISCRIMINATING_RESIDUES,
)
from .geometry import min_distance, virtual_cb

#: Longest plausible peptide-bond length. Author numbering can be contiguous
#: across an unmodelled stretch, so adjacency in numbering is not enough to
#: justify computing a dihedral - the atoms must actually be bonded.
MAX_PEPTIDE_BOND = 2.0

#: Backbone atoms, for separating side-chain atoms from the rest.
_BACKBONE = {"N", "CA", "C", "O", "OXT"}


def phi_psi(
    array: AtomArray, chain_id: str, res_id: int
) -> tuple[float | None, float | None]:
    """Backbone phi/psi in degrees, or ``None`` where undefined.

    ``None`` is returned at chain termini, where a required atom is
    unmodelled, and where the neighbouring residue is numbered consecutively
    but not actually bonded (an unmodelled gap).
    """
    get = lambda c, r, a: structures.atom_coord(array, c, r, a)  # noqa: E731

    n, ca, c = get(chain_id, res_id, "N"), get(chain_id, res_id, "CA"), get(chain_id, res_id, "C")
    if n is None or ca is None or c is None:
        return (None, None)

    phi = None
    prev_c = get(chain_id, res_id - 1, "C")
    if prev_c is not None and np.linalg.norm(prev_c - n) <= MAX_PEPTIDE_BOND:
        phi = float(np.rad2deg(struc.dihedral(prev_c, n, ca, c)))

    psi = None
    next_n = get(chain_id, res_id + 1, "N")
    if next_n is not None and np.linalg.norm(c - next_n) <= MAX_PEPTIDE_BOND:
        psi = float(np.rad2deg(struc.dihedral(n, ca, c, next_n)))

    return (phi, psi)


def peptide_density(array: AtomArray, pdb_id: str, peptide_length: int = 10) -> list[dict]:
    """Per-residue completeness, B-factor and occupancy for every peptide copy.

    Reports the whole-structure mean B-factor alongside, because a residue's
    B only means something relative to the rest of the model.
    """
    structure_mean_b = float(np.mean(array.b_factor))
    rows = []
    for pep in structures.peptide_chains(array):
        for res_id in range(1, peptide_length + 1):
            sel = structures.residue(array, pep, res_id)
            observed = len(sel) > 0
            sidechain = [str(a) for a in sel.atom_name if str(a) not in _BACKBONE]
            rows.append(
                {
                    "pdb_id": pdb_id.upper(),
                    "peptide_chain": pep,
                    "res_id": res_id,
                    "res_name": str(sel.res_name[0]) if observed else "",
                    "observed": observed,
                    "n_atoms": int(len(sel)),
                    "mean_b_factor": round(float(np.mean(sel.b_factor)), 2) if observed else None,
                    "min_occupancy": round(float(np.min(sel.occupancy)), 2) if observed else None,
                    "sidechain_atoms": "|".join(sidechain) if sidechain else "",
                    "structure_mean_b_factor": round(structure_mean_b, 2),
                }
            )
    return rows


def peptide_dihedrals(
    array: AtomArray, pdb_id: str, first: int = 3, last: int = 9
) -> list[dict]:
    """Backbone phi/psi across the central peptide residues, for every copy.

    The default range brackets position 6 with a residue either side, which is
    the region that reorganises between the free and receptor-bound states.
    """
    rows = []
    for pep in structures.peptide_chains(array):
        for res_id in range(first, last + 1):
            sel = structures.residue(array, pep, res_id)
            if len(sel) == 0:
                rows.append(
                    {
                        "pdb_id": pdb_id.upper(),
                        "peptide_chain": pep,
                        "res_id": res_id,
                        "res_name": "",
                        "phi_deg": None,
                        "psi_deg": None,
                        "note": "not modelled",
                    }
                )
                continue
            phi, psi = phi_psi(array, pep, res_id)
            rows.append(
                {
                    "pdb_id": pdb_id.upper(),
                    "peptide_chain": pep,
                    "res_id": res_id,
                    "res_name": str(sel.res_name[0]),
                    "phi_deg": round(phi, 1) if phi is not None else None,
                    "psi_deg": round(psi, 1) if psi is not None else None,
                    "note": (
                        "" if (phi is not None and psi is not None) else "chain break or terminus"
                    ),
                }
            )
    return rows


def sidechain_sasa(
    array: AtomArray,
    chain_id: str,
    res_id: int,
    context_chains: list[str],
) -> float | None:
    """Solvent-accessible area of one side chain, in the context of given chains.

    Context matters: the same residue has a very different exposure with and
    without a receptor bound, and that difference is the whole question here.

    Returns ``None`` for a residue with no side-chain atoms (glycine) or one
    that is not modelled.
    """
    context = array[np.isin(array.chain_id, context_chains)]
    target = (context.chain_id == chain_id) & (context.res_id == res_id)
    sidechain = target & ~np.isin(context.atom_name, list(_BACKBONE))
    if not sidechain.any():
        return None
    per_atom = struc.sasa(context)
    return round(float(np.nansum(per_atom[sidechain])), 2)


def position6_orientation(array: AtomArray, pdb_id: str) -> list[dict]:
    """Which way peptide position 6 points, per copy.

    Distance is measured to the tip of the HLA Arg114 side chain, the partner
    the Asp carboxylate pairs with when buried. Two distances are reported:

    - ``ca_to_ref``: where the backbone sits.
    - ``cb_to_ref``: where the side chain points.

    The informative quantity is their difference. ``cb_to_ref`` greater than
    ``ca_to_ref`` means the side chain leans away from the groove (open);
    smaller means it leans into it (closed).

    For glycine the Cbeta is reconstructed from the backbone frame, which is
    what makes a wild-type structure comparable with a mutant one.
    """
    pairs = structures.pair_peptides_with_heavy_chains(array)
    ref_resname, ref_res_id, ref_atom = GROOVE_REFERENCE_ATOM
    has_tcr = STRUCTURES[pdb_id.upper()].has_tcr

    rows = []
    for pep, heavy in pairs.items():
        sel = structures.residue(array, pep, MUTANT_POSITION)
        if len(sel) == 0:
            rows.append(
                {
                    "pdb_id": pdb_id.upper(),
                    "peptide_chain": pep,
                    "heavy_chain": heavy,
                    "res_name": "",
                    "note": "position 6 not modelled",
                }
            )
            continue

        ref = structures.atom_coord(array, heavy, ref_res_id, ref_atom)
        if ref is None:
            raise ValueError(
                f"{pdb_id} chain {heavy}: {ref_resname}{ref_res_id} {ref_atom} not modelled; "
                "cannot measure groove orientation"
            )

        ca = structures.atom_coord(array, pep, MUTANT_POSITION, "CA")
        real_cb = structures.atom_coord(array, pep, MUTANT_POSITION, "CB")
        n = structures.atom_coord(array, pep, MUTANT_POSITION, "N")
        c = structures.atom_coord(array, pep, MUTANT_POSITION, "C")
        if ca is None or n is None or c is None:
            raise ValueError(f"{pdb_id} chain {pep}: incomplete backbone at position 6")

        cb = real_cb if real_cb is not None else virtual_cb(n, ca, c)
        ca_to_ref = float(np.linalg.norm(ca - ref))
        cb_to_ref = float(np.linalg.norm(cb - ref))

        # SASA in the pMHC context only, so that free and receptor-bound
        # entries are compared on the same basis. The beta-2-microglobulin must
        # be the one belonging to *this* copy; in a multi-copy asymmetric unit
        # several sit within contact distance of the same heavy chain.
        b2m = structures.pair_b2m_with_heavy_chains(array)[heavy]
        context = [pep, heavy, b2m]

        rows.append(
            {
                "pdb_id": pdb_id.upper(),
                "peptide_chain": pep,
                "heavy_chain": heavy,
                "res_name": str(sel.res_name[0]),
                "tcr_present": has_tcr,
                "cb_source": "observed" if real_cb is not None else "reconstructed",
                "ca_to_ref": round(ca_to_ref, 2),
                "cb_to_ref": round(cb_to_ref, 2),
                "cb_minus_ca": round(cb_to_ref - ca_to_ref, 2),
                "points": "out of groove" if cb_to_ref > ca_to_ref else "into groove",
                "sidechain_sasa_pmhc_only": sidechain_sasa(
                    array, pep, MUTANT_POSITION, context
                ),
                "reference": f"{ref_resname}{ref_res_id}:{ref_atom} (chain {heavy})",
                "note": "",
            }
        )
    return rows


def mutant_contacts(array: AtomArray, pdb_id: str) -> list[dict]:
    """Distances from the Asp6 carboxylate to the TCR and HLA residues of interest.

    Requires a modelled Asp at position 6. HLA partners are reported for any
    such entry; TCR partners only where ``constants`` defines them, which is
    7OW6 alone - the per-residue decomposition was done on that complex, and
    its residue numbering does not transfer to a different TCR.

    Raises:
        ValueError: if position 6 is not a modelled aspartate, or if a TCR
            residue does not have the expected identity - the latter would mean
            the alpha/beta chain assignment in ``constants`` is wrong for this
            entry.
    """
    pairs = structures.pair_peptides_with_heavy_chains(array)
    tcr_residues = TCR_DISCRIMINATING_RESIDUES.get(pdb_id.upper(), [])
    rows = []

    for pep, heavy in pairs.items():
        sel = structures.residue(array, pep, MUTANT_POSITION)
        if len(sel) == 0 or str(sel.res_name[0]) != "ASP":
            raise ValueError(
                f"{pdb_id} chain {pep}: position 6 is "
                f"{str(sel.res_name[0]) if len(sel) else 'unmodelled'}, expected ASP"
            )
        carboxylate = sel.coord[np.isin(sel.atom_name, list(ASP_CARBOXYLATE_ATOMS))]
        if len(carboxylate) != len(ASP_CARBOXYLATE_ATOMS):
            raise ValueError(f"{pdb_id} chain {pep}: incomplete Asp6 carboxylate")
        all_atoms = sel.coord

        for chain_id, res_id, expected, label, ddg in tcr_residues:
            partner = structures.residue(array, chain_id, res_id)
            if len(partner) == 0:
                rows.append(
                    {
                        "pdb_id": pdb_id.upper(),
                        "peptide_chain": pep,
                        "partner": label,
                        "partner_chain": chain_id,
                        "partner_res_id": res_id,
                        "partner_res_name": "",
                        "role": "TCR",
                        "ddG_kcal_per_mol": ddg,
                        "note": "not modelled",
                    }
                )
                continue
            observed_name = str(partner.res_name[0])
            if observed_name != expected:
                raise ValueError(
                    f"{pdb_id} chain {chain_id} residue {res_id} is {observed_name}, "
                    f"expected {expected} for {label}: check the alpha/beta chain "
                    "assignment in constants.TCR_DISCRIMINATING_RESIDUES"
                )
            rows.append(
                {
                    "pdb_id": pdb_id.upper(),
                    "peptide_chain": pep,
                    "partner": label,
                    "partner_chain": chain_id,
                    "partner_res_id": res_id,
                    "partner_res_name": observed_name,
                    "role": "TCR",
                    "ddG_kcal_per_mol": ddg,
                    "to_carboxylate": round(min_distance(carboxylate, partner.coord), 2),
                    "to_any_asp6_atom": round(min_distance(all_atoms, partner.coord), 2),
                    "note": "",
                }
            )

        for res_id, expected in sorted(HLA_SURVEY_RESIDUES.items()):
            partner = structures.residue(array, heavy, res_id)
            if len(partner) == 0:
                continue
            observed_name = str(partner.res_name[0])
            rows.append(
                {
                    "pdb_id": pdb_id.upper(),
                    "peptide_chain": pep,
                    "partner": f"{observed_name}{res_id}",
                    "partner_chain": heavy,
                    "partner_res_id": res_id,
                    "partner_res_name": observed_name,
                    "role": "HLA",
                    "ddG_kcal_per_mol": None,
                    "to_carboxylate": round(min_distance(carboxylate, partner.coord), 2),
                    "to_any_asp6_atom": round(min_distance(all_atoms, partner.coord), 2),
                    "note": "" if observed_name == expected else f"expected {expected}",
                }
            )
    return rows
