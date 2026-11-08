"""Tests for the Levy (2010) rSASA interface classification.

This classifier is the project-wide definition of an interface residue, so every
later stage inherits whatever it gets wrong. The structural anchors below are
what make it trustworthy: they are cases where the right answer is known from
independent evidence rather than from this code.
"""

from __future__ import annotations

import biotite.structure as struc
import numpy as np
import pytest

from kras_pmhc import interface, structures
from kras_pmhc.constants import MUTANT_POSITION


def _two_residue_array(separation: float) -> struc.AtomArray:
    """Two single-atom glycines on separate chains, a given distance apart.

    Single-atom residues keep the expected areas analytic: at large separation
    each is a free sphere, and the only thing that can change its area is the
    other atom.
    """
    array = struc.AtomArray(2)
    array.coord = np.array([[0.0, 0.0, 0.0], [separation, 0.0, 0.0]])
    array.chain_id = np.array(["A", "B"])
    array.res_id = np.array([1, 1])
    array.ins_code = np.array(["", ""])
    array.res_name = np.array(["GLY", "GLY"])
    array.atom_name = np.array(["CA", "CA"])
    array.element = np.array(["C", "C"])
    array.hetero = np.array([False, False])
    return array


# --- behaviour at the boundaries ------------------------------------------

def test_chains_far_apart_have_no_interface():
    """Nothing is occluded, so nothing may be classified as interface."""
    records = interface.classify(_two_residue_array(50.0), ["A"], ["B"])

    assert len(records) == 2
    assert not any(r.changed for r in records)
    assert all(r.klass in ("interior", "surface") for r in records)
    assert all(r.delta_rsasa == 0.0 for r in records)
    assert interface.interface_residues(records) == []
    assert interface.interface_area_by_chain(records) == {}


def test_touching_chains_produce_an_interface():
    records = interface.classify(_two_residue_array(3.5), ["A"], ["B"])

    assert all(r.changed for r in records)
    assert len(interface.interface_residues(records)) == 2
    assert all(r.delta_rsasa > 0.0 for r in records)


def test_the_change_tolerance_is_honoured_and_configurable():
    """A residue is interface only if its loss exceeds the tolerance.

    Set the tolerance just above the observed loss and the same structure must
    report no interface - this is the knob, not a fixed ``> 0`` test.
    """
    array = _two_residue_array(3.5)
    loss = max(r.sasa_loss for r in interface.classify(array, ["A"], ["B"]))
    assert loss > interface.DEFAULT_TOLERANCE_A2

    strict = interface.classify(array, ["A"], ["B"], tolerance=loss + 1.0)
    assert not any(r.changed for r in strict)
    assert interface.interface_residues(strict) == []


def test_an_unknown_residue_raises_rather_than_being_skipped():
    array = _two_residue_array(3.5)
    array.res_name = np.array(["UNK", "GLY"])
    with pytest.raises(KeyError, match="UNK"):
        interface.classify(array, ["A"], ["B"])


@pytest.mark.parametrize(
    "group_a,group_b,match",
    [
        (["A"], [], "at least one chain"),
        (["A"], ["A"], "overlap"),
        (["A"], ["Z"], "not in structure"),
    ],
)
def test_group_specification_is_validated(group_a, group_b, match):
    with pytest.raises(ValueError, match=match):
        interface.classify(_two_residue_array(3.5), group_a, group_b)


def test_classes_partition_the_residues():
    """Every residue gets exactly one of the five classes."""
    records = interface.classify(_two_residue_array(3.5), ["A"], ["B"])
    allowed = {"interior", "surface", "support", "rim", "core"}
    assert {r.klass for r in records} <= allowed


# --- structural anchors ---------------------------------------------------

def test_asp6_is_core_against_the_hla_groove(ow6):
    """The known-correct answer, from independent evidence.

    Stage 01 measured the Asp6 side chain at 7.4 A^2 in this structure against
    53-82 A^2 on the open backbone, so the groove demonstrably buries it. A
    classifier that does not call this ``core`` is wrong.
    """
    peptide = structures.peptide_chains(ow6)
    records = interface.classify(ow6, peptide, structures.heavy_chains(ow6))

    asp6 = interface.find(records, peptide[0], MUTANT_POSITION)
    assert asp6 is not None
    assert asp6.res_name == "ASP"
    assert asp6.klass == "core"
    assert asp6.rsasa_monomer >= interface.RSASA_THRESHOLD
    assert asp6.rsasa_complex < interface.RSASA_THRESHOLD
    assert asp6 in interface.buried_interface_residues(records)


def test_asp6_is_not_buried_by_the_tcr(ow6):
    """The TCR reads the aspartate without desolvating it.

    Poole et al. attribute discrimination to long-range electrostatics rather
    than contacts, and stage 01 measured the dominant energetic contributors at
    13-16 A from the carboxylate. So the aspartate must *not* come out buried
    against the TCR - and this is what the two-group parameter exists for: the
    per-chain reading would have the HLA's burial contaminate this answer.
    """
    peptide = structures.peptide_chains(ow6)
    tcr = ["D", "E"]
    records = interface.classify(ow6, peptide, tcr)

    asp6 = interface.find(records, peptide[0], MUTANT_POSITION)
    assert asp6 is not None
    assert asp6.klass not in interface.BURIED_INTERFACE_CLASSES
    assert asp6.rsasa_complex >= interface.RSASA_THRESHOLD


def test_grouping_changes_the_answer_for_the_same_residue(ow6):
    """Documents why ``classify`` takes groups rather than single chains."""
    peptide = structures.peptide_chains(ow6)
    position = (peptide[0], MUTANT_POSITION)

    against_hla = interface.find(
        interface.classify(ow6, peptide, structures.heavy_chains(ow6)), *position
    )
    against_tcr = interface.find(
        interface.classify(ow6, peptide, ["D", "E"]), *position
    )

    assert against_hla.rsasa_monomer == against_tcr.rsasa_monomer
    assert against_hla.sasa_loss > 4 * against_tcr.sasa_loss


def test_only_the_named_chains_take_part(ow6):
    """Chains outside both groups must not affect either calculation.

    Otherwise an interface would silently depend on what else happened to be in
    the file, and the same two chains would classify differently in the full
    complex than in an extracted pair.
    """
    peptide = structures.peptide_chains(ow6)
    heavy = structures.heavy_chains(ow6)

    full = interface.classify(ow6, peptide, heavy)
    extracted = interface.classify(
        ow6[np.isin(ow6.chain_id, peptide + heavy)], peptide, heavy
    )

    assert [r.klass for r in full] == [r.klass for r in extracted]


def test_interface_area_by_chain_splits_a_multichain_partner(ow6):
    """The metric the stage-05 pilot gate reads.

    With the peptide on one side and the two TCR chains on the other, the
    per-chain totals must account for all three chains and be positive.
    """
    peptide = structures.peptide_chains(ow6)
    records = interface.classify(ow6, peptide, ["D", "E"])
    areas = interface.interface_area_by_chain(records)

    assert set(areas) == {peptide[0], "D", "E"}
    assert all(area > 0.0 for area in areas.values())
