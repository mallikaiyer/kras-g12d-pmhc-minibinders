"""Tests that encode the structural claims the design plan rests on.

These are not unit tests of convenience helpers. Each one asserts a published
or measured fact about the target, so that if a refactor, a biotite upgrade or
a re-download changes the numbers, the suite says so rather than letting a
revised figure propagate quietly into the methods document.
"""

from __future__ import annotations

import pytest

from kras_pmhc import structures, survey
from kras_pmhc.constants import MUTANT_POSITION

# --- chain interpretation -------------------------------------------------

def test_peptide_identity_holds_for_every_entry(ow3, ow4, ow5, ow6, pb2):
    """Every modelled peptide position matches the expected sequence.

    This is the guard against analysing the wrong chain: it would fail loudly
    if chain-role assignment picked up, say, a TCR CDR loop as the peptide.
    """
    for pdb_id, array in [
        ("7OW3", ow3), ("7OW4", ow4), ("7OW5", ow5), ("7OW6", ow6), ("7PB2", pb2)
    ]:
        observed = structures.check_peptide_identity(array, pdb_id)
        assert observed, f"{pdb_id}: no peptide chain found"


def test_peptide_pairing_is_derived_not_assumed(ow3):
    """7OW3 has four copies; each peptide pairs with its own heavy chain."""
    pairs = structures.pair_peptides_with_heavy_chains(ow3)
    assert pairs == {"C": "A", "F": "D", "I": "G", "L": "J"}


# --- completeness ---------------------------------------------------------

def test_7ow4_has_no_position6_in_any_copy(ow4):
    """The G12D free structure never resolves the mutated residue.

    This is why Route A cannot simply use the mutant's own free structure, and
    why the open backbone has to be taken from the wild-type entry instead.
    """
    rows = survey.peptide_density(ow4, "7OW4")
    position6 = [r for r in rows if r["res_id"] == MUTANT_POSITION]
    assert len(position6) == 4, "expected four copies in the asymmetric unit"
    assert not any(r["observed"] for r in position6)


def test_7ow3_has_exactly_two_complete_copies(ow3):
    """Two of four wild-type copies resolve the central residues 5-7.

    Route A's open backbone comes from these. The other two are gapped, which
    is what the paper means by density being "only partially observed".
    """
    rows = survey.peptide_density(ow3, "7OW3")
    complete = {
        chain
        for chain in structures.peptide_chains(ow3)
        if all(
            r["observed"]
            for r in rows
            if r["peptide_chain"] == chain and r["res_id"] in (5, 6, 7)
        )
    }
    assert complete == {"C", "L"}


# --- the Route A gate -----------------------------------------------------

def test_open_backbone_is_not_glycine_only(ow3):
    """Position 6 in the free wild-type structure has negative phi.

    Glycine-exclusive Ramachandran territory is overwhelmingly positive phi.
    A positive value here would mean an aspartate could not occupy the open
    backbone, which would eliminate Route A outright.
    """
    rows = survey.peptide_dihedrals(ow3, "7OW3")
    position6 = [
        r for r in rows
        if r["res_id"] == MUTANT_POSITION and r["phi_deg"] is not None
    ]
    assert len(position6) == 2, "expected phi for the two complete copies"
    for row in position6:
        assert row["phi_deg"] < 0, (
            f"chain {row['peptide_chain']} phi={row['phi_deg']} is positive: "
            "glycine-only backbone would invalidate Route A"
        )


# --- the open-versus-closed measurement -----------------------------------

def _orientation_by_chain(array, pdb_id):
    return {
        r["peptide_chain"]: r
        for r in survey.position6_orientation(array, pdb_id)
        if not r.get("note")
    }


def test_position6_points_out_of_the_groove_when_free(ow3):
    rows = _orientation_by_chain(ow3, "7OW3")
    assert set(rows) == {"C", "L"}
    for chain, row in rows.items():
        assert row["cb_minus_ca"] > 0.5, f"chain {chain} should lean away from the groove"
        assert row["points"] == "out of groove"


@pytest.mark.parametrize("fixture_name,pdb_id", [("ow5", "7OW5"), ("ow6", "7OW6")])
def test_position6_points_into_the_groove_when_receptor_bound(
    request, fixture_name, pdb_id
):
    """Both the wild-type and the mutant close on TCR binding.

    7OW5 is the important one: it is the same glycine as 7OW3 and the same
    reconstruction method, so the comparison between them isolates the
    conformational change from the identity of the residue.
    """
    array = request.getfixturevalue(fixture_name)
    rows = _orientation_by_chain(array, pdb_id)
    for chain, row in rows.items():
        assert row["cb_minus_ca"] < -0.5, f"{pdb_id} chain {chain} should lean into the groove"
        assert row["points"] == "into groove"


def test_the_backbone_itself_moves_between_states(ow3, ow6):
    """The Ca of position 6 sits ~4 A further from the groove when free.

    Guards against the orientation flip being an artefact of the Cbeta
    reconstruction alone: the backbone moves too.
    """
    free = _orientation_by_chain(ow3, "7OW3")
    bound = _orientation_by_chain(ow6, "7OW6")
    free_mean = sum(r["ca_to_ref"] for r in free.values()) / len(free)
    bound_mean = sum(r["ca_to_ref"] for r in bound.values()) / len(bound)
    assert free_mean - bound_mean > 3.5, (
        f"free Ca {free_mean:.2f} A vs bound {bound_mean:.2f} A: "
        "expected a substantial backbone shift"
    )


# --- the mechanism claim --------------------------------------------------

def test_asp6_is_engaged_by_the_hla_in_the_closed_state(ow6):
    """The carboxylate salt-bridges Arg114 and hydrogen-bonds Gln70.

    These are the interactions a designed binder would have to outcompete,
    and they are why the closed state is not a viable design target.
    """
    rows = {r["partner"]: r for r in survey.mutant_contacts(ow6, "7OW6")}
    assert rows["ARG114"]["to_carboxylate"] < 3.2
    assert rows["GLN70"]["to_carboxylate"] < 3.2


def test_the_discriminating_tcr_residues_are_distant(ow6):
    """The two largest energetic contributors are >13 A from the carboxylate.

    This is the evidence that the TCR's specificity is long-range electrostatic
    rather than contact-mediated, and therefore that Route B cannot be designed
    by placing contacts.
    """
    rows = {r["partner"]: r for r in survey.mutant_contacts(ow6, "7OW6")}
    assert rows["R50a"]["to_carboxylate"] > 13.0
    assert rows["K70a"]["to_carboxylate"] > 13.0
    # Only one named residue comes near, and only via a non-carboxylate atom.
    assert rows["H100b"]["to_carboxylate"] > 6.0


def test_tcr_table_is_not_applied_to_a_different_tcr(pb2):
    """7PB2 carries the parent TCR, so no TCR rows should be produced.

    The residue numbering in the decomposition table belongs to JDIa41b1.
    Applying it here would report distances for unrelated residues; this test
    pins the scoping that prevents that.
    """
    rows = survey.mutant_contacts(pb2, "7PB2")
    assert rows, "HLA contacts should still be reported"
    assert all(r["role"] == "HLA" for r in rows)
