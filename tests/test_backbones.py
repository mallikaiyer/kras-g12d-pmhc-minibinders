"""Tests for generated-backbone metrics and the stage-05 manifest.

The synthetic binder below is not a design and its numbers mean nothing
biologically. Its purpose is to exercise every code path the ingest script takes
*before* GPU time is spent, because discovering a parsing bug after a six-hour
generation run would be expensive. It is a real 70-residue protein fragment
translated into contact with the peptide, so chain handling, residue numbering,
the Levy classification and the secondary-structure annotation all see
well-formed input.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from biotite.structure.io.pdb import PDBFile

from kras_pmhc import backbones, interface
from kras_pmhc.constants import ASP_CARBOXYLATE_ATOMS, MUTANT_POSITION

ROOT = Path(__file__).resolve().parents[1]
TARGET = (
    ROOT / "data" / "targets_relaxed"
    / "7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb"
)
MANIFEST = ROOT / "data" / "backbones" / "MANIFEST.tsv"

TARGET_CHAINS = ["A", "C"]
FRAGMENT_RANGE = (60, 129)  # 70 residues of chain A, reused as a pseudo-binder


@pytest.fixture(scope="module")
def target():
    if not TARGET.exists():
        pytest.skip("stage 04 target not built; run scripts/04_relax_targets.py")
    return PDBFile.read(str(TARGET)).get_structure(model=1)


#: Contact separation the fixture docks to, A. Close enough to bury surface,
#: far enough to be a contact rather than an overlap.
CONTACT_DISTANCE = 3.5


def _fragment(target, low, high):
    """A chain-A slice reused as a pseudo-binder: real geometry, new chain."""
    mask = (
        (target.chain_id == "A")
        & (target.res_id >= low)
        & (target.res_id <= high)
    )
    fragment = target[mask].copy()
    fragment.chain_id = np.full(fragment.array_length(), "B")
    fragment.res_id = fragment.res_id - low + 1
    return fragment


def _dock(target, fragment, separation=CONTACT_DISTANCE):
    """Slide the fragment outward until it just touches the target.

    Bisection is on the distance to the **whole target**, not to Asp6. Bisecting
    on the carboxylate alone lands the fragment inside the HLA body - a
    comfortable distance from the aspartate while interpenetrating everything
    else - which is a real failure mode and no use as a positive control.
    """
    target_coords = target.coord
    origin = fragment.coord.mean(axis=0)
    direction = target_coords.mean(axis=0) - origin
    direction /= np.linalg.norm(direction)

    def closest(shift):
        moved = fragment.coord + direction * shift
        return np.sqrt(
            ((moved[:, None, :] - target_coords[None, :, :]) ** 2).sum(-1)
        ).min()

    # Walk outward (negative shift) until clear, then bisect back to contact.
    far = 0.0
    while closest(far) < separation:
        far -= 5.0
    lo, hi = far, 0.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if closest(mid) >= separation:
            lo = mid
        else:
            hi = mid

    docked = fragment.copy()
    docked.coord = fragment.coord + direction * lo
    return docked


@pytest.fixture(scope="module")
def synthetic(target):
    """Target plus a 70-residue chain B docked into surface contact with it."""
    return target + _dock(target, _fragment(target, *FRAGMENT_RANGE))


# --- binder_chain ---------------------------------------------------------

def test_binder_chain_found_by_exclusion(synthetic):
    assert backbones.binder_chain(synthetic, TARGET_CHAINS) == "B"


def test_binder_chain_rejects_a_missing_target_chain(synthetic):
    with pytest.raises(ValueError, match="absent"):
        backbones.binder_chain(synthetic, ["A", "C", "Z"])


def test_binder_chain_rejects_an_ambiguous_result(synthetic):
    """No generated chain, and more than one, are both failures."""
    with pytest.raises(ValueError, match="exactly one generated chain"):
        backbones.binder_chain(synthetic, ["A", "B", "C"])
    with pytest.raises(ValueError, match="exactly one generated chain"):
        backbones.binder_chain(synthetic, ["A"])


# --- target_backbone_rmsd -------------------------------------------------

def test_an_untouched_target_has_zero_rmsd(synthetic, target):
    """The check that catches a misparse or a wrong contig.

    Computed without superposition, so a rigid-body shift of the target would
    show up here rather than being absorbed by a fit.
    """
    rmsd = backbones.target_backbone_rmsd(synthetic, target, TARGET_CHAINS)
    assert rmsd == pytest.approx(0.0, abs=1e-6)


def test_a_shifted_target_is_detected(synthetic, target):
    moved = synthetic.copy()
    moved.coord[np.isin(moved.chain_id, TARGET_CHAINS)] += np.array([1.0, 0.0, 0.0])
    assert backbones.target_backbone_rmsd(
        moved, target, TARGET_CHAINS
    ) == pytest.approx(1.0, abs=1e-6)


def test_rmsd_raises_when_nothing_corresponds(synthetic, target):
    renamed = target.copy()
    renamed.chain_id = np.full(renamed.array_length(), "Z")
    with pytest.raises(ValueError, match="no target backbone atoms"):
        backbones.target_backbone_rmsd(synthetic, renamed, TARGET_CHAINS)


# --- shape metrics --------------------------------------------------------

def test_chain_length_counts_residues(synthetic):
    low, high = FRAGMENT_RANGE
    assert backbones.chain_length(synthetic, "B") == high - low + 1
    with pytest.raises(ValueError, match="no atoms"):
        backbones.chain_length(synthetic, "Q")


def test_secondary_structure_fractions_sum_to_one(synthetic):
    sse = backbones.secondary_structure_fraction(synthetic, "B")
    assert set(sse) == {"helix", "strand", "coil"}
    assert sum(sse.values()) == pytest.approx(1.0, abs=1e-3)
    assert all(0.0 <= value <= 1.0 for value in sse.values())


# --- the ingest path ------------------------------------------------------

def test_the_docked_binder_does_not_overlap_the_target(synthetic):
    """The fixture must be a contact, not an interpenetration."""
    binder = synthetic.coord[synthetic.chain_id == "B"]
    target_coords = synthetic.coord[np.isin(synthetic.chain_id, TARGET_CHAINS)]
    separation = np.sqrt(
        ((binder[:, None, :] - target_coords[None, :, :]) ** 2).sum(-1)
    ).min()
    assert separation == pytest.approx(CONTACT_DISTANCE, abs=0.1)


def test_the_synthetic_binder_produces_a_usable_interface(synthetic):
    """The exact call stage 05's ingest makes, on known-well-formed input.

    Asserts the classification is non-degenerate and internally consistent: a
    contacting 70-residue chain must yield interface residues on both sides, and
    the per-chain areas must cover the chains in contact.
    """
    records = interface.classify(synthetic, ["B"], TARGET_CHAINS)
    found = interface.interface_residues(records)
    assert found, "a chain in contact must produce an interface"

    sides = {r.group for r in found}
    assert sides == {"a", "b"}, "both sides of an interface must be populated"

    areas = interface.interface_area_by_chain(records)
    assert "B" in areas and areas["B"] > 0.0
    assert set(areas) <= {"A", "B", "C"}

    asp6 = interface.find(records, "C", MUTANT_POSITION)
    assert asp6 is not None and asp6.res_name == "ASP"


def test_an_interpenetrating_chain_gives_a_one_sided_interface(target):
    """The invariant the ingest script validates, and why it is needed.

    Constructed the naive way - bisecting on distance to the Asp6 carboxylate
    rather than to the whole target - which is how the bug was originally
    found: it parks the fragment 8 A from the aspartate and simultaneously
    inside the HLA body. There it loses area to cavity-lining atoms that had
    almost none to begin with, so the binder side shows an interface and the
    target side shows none. Occlusion is mutual, so that asymmetry can only
    mean interpenetration, which is exactly what the check detects.
    """
    fragment = _fragment(target, 20, 84)
    carboxylate = target.coord[
        (target.chain_id == "C")
        & (target.res_id == MUTANT_POSITION)
        & np.isin(target.atom_name, list(ASP_CARBOXYLATE_ATOMS))
    ]
    origin = fragment.coord.mean(axis=0)
    direction = carboxylate.mean(axis=0) - origin
    direction /= np.linalg.norm(direction)

    def closest(shift):
        moved = fragment.coord + direction * shift
        return np.sqrt(
            ((moved[:, None, :] - carboxylate[None, :, :]) ** 2).sum(-1)
        ).min()

    lo, hi = 0.0, float(np.linalg.norm(carboxylate.mean(axis=0) - origin)) + 20.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if closest(mid) > 8.0 else (lo, mid)
    buried = fragment.copy()
    buried.coord = fragment.coord + direction * lo

    records = interface.classify(target + buried, ["B"], TARGET_CHAINS)
    binder_side = [
        r for r in interface.interface_residues(records) if r.chain_id == "B"
    ]
    target_side = [
        r for r in interface.interface_residues(records) if r.chain_id != "B"
    ]
    assert binder_side and not target_side


# --- outputs, once stage 05 has run --------------------------------------

def _manifest_rows():
    if not MANIFEST.exists():
        pytest.skip("stage 05 not ingested; run scripts/05_ingest_backbones.py")
    with MANIFEST.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_every_backbone_left_the_target_alone():
    """Mirrors the script's threshold deliberately: if the script loosens its
    bound, this should fail rather than follow."""
    for row in _manifest_rows():
        rmsd = float(row["target_rmsd"])
        assert rmsd < 0.1, f"{row['backbone_id']}: target moved {rmsd} A"


def test_every_binder_is_the_requested_length():
    for row in _manifest_rows():
        length = int(row["binder_length"])
        assert 60 <= length <= 75, f"{row['backbone_id']}: {length} residues"


def test_no_binder_overlaps_the_target():
    for row in _manifest_rows():
        closest = float(row["min_dist_to_target"])
        assert closest >= 2.0, f"{row['backbone_id']}: closest {closest} A"


def test_every_interface_is_mutual():
    """One-sided occlusion is impossible; seeing it means interpenetration."""
    for row in _manifest_rows():
        binder = int(row["n_interface_binder"])
        target_side = int(row["n_interface_target"])
        assert bool(binder) == bool(target_side), (
            f"{row['backbone_id']}: one-sided interface "
            f"(binder {binder}, target {target_side})"
        )


def test_the_gate_column_agrees_with_its_two_criteria():
    """The recorded verdict must be reproducible from the recorded evidence."""
    for row in _manifest_rows():
        expected = (
            row["asp6_class"] in interface.INTERFACE_CLASSES
            and float(row["peptide_share"]) >= 0.20
        )
        assert (row["passes_gate"] == "True") is expected, row["backbone_id"]
