"""Tests for protonation and restrained minimisation.

The restraint is what preserves Route A's central claim - that the peptide
backbone is an experimental observation rather than a model - so what it acts
on is worth pinning. The manifest assertions guard the relaxation outputs
against silently drifting if the protocol or the force field changes.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from kras_pmhc import relax

TARGETS = Path(__file__).resolve().parents[1] / "data" / "targets"
RELAXED = Path(__file__).resolve().parents[1] / "data" / "targets_relaxed"

# Thresholds mirror scripts/04_relax_targets.py. Duplicated deliberately: if
# the script loosens a bound, these tests should fail rather than follow.
MAX_BACKBONE_RMSD = 0.3
MAX_CHI1_DRIFT = 30.0
MIN_CARBOXYLATE_SASA = 25.0
MIN_HEAVY_ATOM_APPROACH = 2.4


@pytest.fixture(scope="module")
def protonated():
    path = TARGETS / "7ow3_L_G6D_chi1m70_chi2p120.pdb"
    if not path.exists():
        pytest.skip("stage 03 targets not built; run scripts/03_build_targets.py")
    return relax.protonate(path)


def test_protonation_adds_hydrogens(protonated):
    topology, positions = protonated
    elements = [a.element.symbol for a in topology.atoms() if a.element is not None]
    assert elements.count("H") > 0
    assert topology.getNumAtoms() == len(positions)


def test_protonation_builds_no_missing_residues(protonated):
    """The target must contain no invented backbone.

    Copy C of 7OW3 starts at heavy-chain residue 2. PDBFixer would happily
    build residue 1, and `protonate` deliberately prevents that.
    """
    topology, _ = protonated
    # The deposited copy has 385 residues; protonation adds atoms, never residues.
    assert topology.getNumResidues() == 385


def test_restraint_acts_on_backbone_heavy_atoms_only(protonated):
    """Exactly N, CA, C and O - no hydrogens, no side chains, no OXT.

    OXT is excluded because a terminal carboxylate oxygen's placement is a
    modelling choice rather than part of the observed trace.
    """
    topology, positions = protonated
    _, restrained = relax.backbone_restraint_force(topology, positions)

    atoms = list(topology.atoms())
    restrained_names = {atoms[i].name for i in restrained}
    assert restrained_names == set(relax.BACKBONE_ATOMS)

    expected = sum(1 for a in atoms if a.name in relax.BACKBONE_ATOMS)
    assert len(restrained) == expected
    assert "OXT" not in restrained_names


# --- outputs --------------------------------------------------------------

def _manifest_rows():
    path = RELAXED / "MANIFEST.tsv"
    if not path.exists():
        pytest.skip("stage 04 not run; run scripts/04_relax_targets.py")
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_every_structure_kept_its_backbone():
    for row in _manifest_rows():
        rmsd = float(row["backbone_rmsd"])
        assert rmsd < MAX_BACKBONE_RMSD, f"{row['target_id']}: backbone RMSD {rmsd}"


def test_every_mutant_stayed_in_its_rotamer_well():
    rows = [r for r in _manifest_rows() if r["is_mutant"] == "True"]
    assert len(rows) == 5, f"expected five mutants, found {len(rows)}"
    for row in rows:
        before, after = float(row["chi1_before"]), float(row["chi1_after"])
        drift = abs((before - after + 180.0) % 360.0 - 180.0)
        assert drift < MAX_CHI1_DRIFT, f"{row['target_id']}: chi1 moved {drift:.1f} deg"


def test_every_mutant_kept_an_exposed_carboxylate():
    """Relaxation must not bury the group a binder is meant to engage."""
    for row in [r for r in _manifest_rows() if r["is_mutant"] == "True"]:
        sasa = float(row["carboxylate_sasa"])
        assert sasa >= MIN_CARBOXYLATE_SASA, f"{row['target_id']}: SASA {sasa} A^2"


def test_every_mutant_is_sterically_clean():
    for row in [r for r in _manifest_rows() if r["is_mutant"] == "True"]:
        closest = float(row["closest_heavy_atom"])
        assert closest >= MIN_HEAVY_ATOM_APPROACH, (
            f"{row['target_id']}: closest approach {closest} A"
        )


def test_phi_stayed_negative_for_every_peptide_bearing_structure():
    """The stage-02 Ramachandran gate must survive relaxation.

    Glycine-exclusive backbone is overwhelmingly positive phi; a relaxed
    positive phi would mean the open backbone no longer admits an aspartate.
    """
    for row in _manifest_rows():
        if not row["phi_after"]:
            continue
        assert float(row["phi_after"]) < 0.0, f"{row['target_id']}: phi became positive"


def test_the_groove_does_not_reach_for_the_aspartate():
    """Arg114 must relax by the same amount with and without the aspartate.

    This is the one piece of physical evidence stage 04 contributes. If the
    groove pulled measurably toward an open Asp6, that would be a warning for
    the conformational-selection premise. The wild-type and apo references are
    what make the mutant number interpretable, which is why all nine
    structures get an identical protocol.
    """
    rows = _manifest_rows()
    mutants = [float(r["sidechain_rmsd_ARG114"]) for r in rows if r["is_mutant"] == "True"]
    controls = [float(r["sidechain_rmsd_ARG114"]) for r in rows if r["is_mutant"] != "True"]
    assert mutants and controls

    # Mutant displacement must not exceed the controls by a meaningful margin.
    assert max(mutants) - max(controls) < 0.15, (
        f"Arg114 moves more with the aspartate present "
        f"(mutants up to {max(mutants):.3f} A, controls up to {max(controls):.3f} A)"
    )
