"""Checks on the literature and challenge constants.

The point of these is that a constant copied from a web page is the one kind of
input nothing else in the pipeline can validate. The heavy-chain sequence can
be validated, though, because the design target is a crop of it - so it is,
here, against the structure every other stage was built from.
"""

from __future__ import annotations

import biotite.structure as struc
import pytest
from biotite.sequence import ProteinSequence

from kras_pmhc import constants

TARGET = "7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb"


def test_the_heavy_chain_is_the_stated_length():
    assert len(constants.HLA_A1101_HEAVY_CHAIN) == 281


def test_the_tien_table_covers_the_standard_twenty():
    assert len(constants.MAX_ASA_TIEN_2013) == 20
    assert constants.MAX_ASA_TIEN_2013["GLY"] < constants.MAX_ASA_TIEN_2013["TRP"]


@pytest.mark.needs_data
def test_the_heavy_chain_matches_the_design_target_exactly(relaxed_target):
    """Confirms the allele and, more usefully, the residue numbering.

    Every residue index in this project - the Asp6 hotspot, the A1-180 contig,
    the Levy tables - is expressed in the design target's numbering. This test
    is what ties that numbering to the sequence the challenge specifies.
    """
    chain_a = relaxed_target[relaxed_target.chain_id == "A"]
    starts = struc.get_residue_starts(chain_a)
    observed = "".join(
        ProteinSequence.convert_letter_3to1(name) for name in chain_a.res_name[starts]
    )

    assert len(observed) == constants.HLA_DESIGN_TARGET_RESIDUES
    expected = constants.HLA_A1101_HEAVY_CHAIN[: len(observed)]
    assert observed == expected, "the design target is not a prefix of the heavy chain"


@pytest.mark.needs_data
def test_the_peptide_is_the_g12d_decamer(relaxed_target):
    peptide = relaxed_target[relaxed_target.chain_id == "C"]
    starts = struc.get_residue_starts(peptide)
    observed = "".join(
        ProteinSequence.convert_letter_3to1(name) for name in peptide.res_name[starts]
    )
    assert observed == "VVVGADGVGK"


def test_beta_2_microglobulin_is_the_stated_length():
    assert len(constants.B2M_SEQUENCE) == 99


@pytest.mark.needs_data
def test_beta_2_microglobulin_matches_the_crystal_structure(ow3):
    """Residues 1-99 of 7OW3 chain B, past the initiator Met at residue 0.

    The offset is the point of this test: the challenge gives the mature chain
    and the crystal an expressed one, and stage 11 must describe the construct
    the challenge asked about rather than the one that was crystallised.
    """
    chain_b = ow3[(ow3.chain_id == "B") & struc.filter_amino_acids(ow3)]
    starts = struc.get_residue_starts(chain_b)
    observed = "".join(
        ProteinSequence.convert_letter_3to1(name) for name in chain_b.res_name[starts]
    )
    res_ids = chain_b.res_id[starts]

    assert res_ids.min() == 0, "expected the initiator Met at residue 0"
    assert observed[0] == "M"
    assert constants.B2M_STRUCTURE_HAS_INITIATOR_MET is True
    assert observed[1:] == constants.B2M_SEQUENCE
