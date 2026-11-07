"""Tests for biological-assembly extraction and the Gly6->Asp mutation.

The assembly tests matter because beta-2-microglobulin cannot be assigned to a
heavy chain by proximity in a multi-copy asymmetric unit, and getting it wrong
would mean simulating a complex assembled from two different copies.

The mutation tests matter because Route A's whole claim is that the backbone is
an experimental observation rather than a model. If the mutation perturbed the
backbone, that claim would no longer hold.
"""

from __future__ import annotations

import biotite.structure as struc
import numpy as np
import pytest

from kras_pmhc import structures
from kras_pmhc.build import ASP_ADDED_ATOMS, mutate_gly_to_asp
from kras_pmhc.constants import MUTANT_POSITION

BACKBONE = ("N", "CA", "C", "O")


# --- biological assembly --------------------------------------------------

def test_b2m_pairing_on_a_four_copy_asymmetric_unit(ow3):
    """7OW3's four copies pair A-B, D-E, G-H, J-K.

    Closest approach cannot produce this: chains B, E and K all sit within
    3.3 A of heavy chain J. Contact count separates them by roughly tenfold.
    """
    assert structures.pair_b2m_with_heavy_chains(ow3) == {
        "A": "B", "D": "E", "G": "H", "J": "K"
    }


def test_b2m_pairing_on_a_two_copy_asymmetric_unit(pb2):
    """7PB2 pairs A-B and F-G. Taking the first candidate would give F-B."""
    assert structures.pair_b2m_with_heavy_chains(pb2) == {"A": "B", "F": "G"}


def test_pairing_refuses_a_structure_with_no_b2m(ow3):
    """Absent a partner the function raises rather than inventing one."""
    heavy_only = ow3[~np.isin(ow3.chain_id, structures.b2m_chains(ow3))]
    with pytest.raises(ValueError, match="no beta-2-microglobulin"):
        structures.pair_b2m_with_heavy_chains(heavy_only)


def test_extract_complex_returns_one_biological_pmhc(ow3):
    """Peptide L belongs with heavy chain J and beta-2-microglobulin K."""
    complex_ = structures.extract_complex(ow3, "L")
    assert set(np.unique(complex_.chain_id)) == {"L", "J", "K"}

    counts = structures.chain_residue_counts(complex_)
    assert counts["L"] == 10          # the peptide
    assert counts["K"] == 100         # beta-2-microglobulin
    assert counts["J"] >= 250         # class I heavy chain


def test_extract_complex_rejects_a_non_peptide_chain(ow3):
    with pytest.raises(ValueError, match="not a peptide chain"):
        structures.extract_complex(ow3, "A")


def test_select_atoms_raises_when_nothing_matches(ow3):
    """Position 6 of a wild-type copy is a glycine, so it has no carboxylate."""
    with pytest.raises(ValueError, match="none of the atoms"):
        structures.select_atoms(ow3, "C", MUTANT_POSITION, ("OD1", "OD2"))


# --- the mutation ---------------------------------------------------------

@pytest.fixture(scope="module")
def mutated(ow3):
    """Copy L mutated at the most populated, most exposed rotamer."""
    complex_ = structures.extract_complex(ow3, "L")
    return complex_, mutate_gly_to_asp(complex_, "L", MUTANT_POSITION, -70.0, 120.0)


def test_mutation_adds_exactly_four_atoms(mutated):
    before, after = mutated
    assert len(after) == len(before) + len(ASP_ADDED_ATOMS)

    added = structures.residue(after, "L", MUTANT_POSITION)
    assert str(added.res_name[0]) == "ASP"
    assert set(ASP_ADDED_ATOMS) <= {str(name) for name in added.atom_name}


def test_mutation_leaves_every_original_coordinate_untouched(mutated):
    """Not merely the mutated residue: nothing observed may move."""
    before, after = mutated
    key = lambda a, i: (str(a.chain_id[i]), int(a.res_id[i]), str(a.atom_name[i]))  # noqa: E731

    after_coords = {key(after, i): after.coord[i] for i in range(len(after))}
    for i in range(len(before)):
        k = key(before, i)
        assert k in after_coords, f"{k} disappeared"
        assert np.array_equal(before.coord[i], after_coords[k]), f"{k} moved"


def test_mutation_preserves_the_backbone_dihedrals(mutated):
    """phi/psi at position 6 must be bit-identical, not merely close."""
    before, after = mutated
    from kras_pmhc.survey import phi_psi

    assert phi_psi(before, "L", MUTANT_POSITION) == phi_psi(after, "L", MUTANT_POSITION)


def test_mutation_sets_the_requested_torsions(mutated):
    _, after = mutated
    get = lambda name: structures.atom_coord(after, "L", MUTANT_POSITION, name)  # noqa: E731

    chi1 = float(np.rad2deg(struc.dihedral(get("N"), get("CA"), get("CB"), get("CG"))))
    chi2 = float(np.rad2deg(struc.dihedral(get("CA"), get("CB"), get("CG"), get("OD1"))))
    assert chi1 == pytest.approx(-70.0, abs=0.1)
    assert chi2 == pytest.approx(120.0, abs=0.1)


def test_mutation_rejects_a_residue_that_is_not_glycine(ow6):
    """Position 6 of 7OW6 is already an aspartate."""
    with pytest.raises(ValueError, match="expected GLY"):
        mutate_gly_to_asp(ow6, "C", MUTANT_POSITION, -70.0, 120.0)


def test_mutation_rejects_an_unmodelled_residue(ow3):
    """Copy F has no coordinates at position 6."""
    with pytest.raises(ValueError, match="not modelled"):
        mutate_gly_to_asp(ow3, "F", MUTANT_POSITION, -70.0, 120.0)
