"""Metrics for generated binder backbones.

RFdiffusion returns backbone-only coordinates (N, CA, C, O) for a new chain
alongside the target it was conditioned on. This module answers the three
questions that have to be settled before a backbone is worth designing a
sequence onto: is it the chain we think it is, did the generator leave the
target alone, and what fold did it produce.

All interface work is delegated to :mod:`kras_pmhc.interface`, which implements
the Levy (2010) rSASA classification used project-wide. Nothing here measures an
interface by distance.
"""

from __future__ import annotations

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray

#: Backbone heavy atoms RFdiffusion emits, and the only ones present on a
#: generated chain. Matching on these makes the target comparison insensitive
#: to whether side chains survived the round trip.
BACKBONE_ATOMS = ("N", "CA", "C", "O")

#: The atoms that define a residue's rigid frame. RFdiffusion diffuses frames
#: and rebuilds the carbonyl oxygen from idealised internal geometry, so O in
#: its output is *not* the input oxygen: measured against the stage-04 target on
#: 9 Nov 2026, N/CA/C agree to ~0.1 A rms while O differs by 0.37 A rms and up
#: to 2.4 A. Comparisons of target geometry therefore use these three.
FRAME_ATOMS = ("N", "CA", "C")

#: Secondary-structure codes from ``biotite.structure.annotate_sse``.
SSE_HELIX = "a"
SSE_STRAND = "b"


def binder_chain(array: AtomArray, target_chains: list[str]) -> str:
    """The generated chain, identified by excluding the target chains.

    Identification by exclusion rather than by an expected letter: RFdiffusion
    and the ColabDesign post-processing both assign chain ids, and hard-coding
    one would fail silently if either changed.

    Raises:
        ValueError: if no chain or more than one chain remains, since both mean
            the output does not have the shape this stage assumes.
    """
    present = set(np.unique(array.chain_id).tolist())
    missing = set(target_chains) - present
    if missing:
        raise ValueError(
            f"target chain(s) {sorted(missing)} absent from the generated "
            f"structure; present: {sorted(present)}"
        )

    remaining = sorted(present - set(target_chains))
    if len(remaining) != 1:
        raise ValueError(
            f"expected exactly one generated chain, found {remaining} "
            f"after excluding target chains {sorted(target_chains)}"
        )
    return remaining[0]


def _backbone_index(
    array: AtomArray, chains: list[str], atoms: tuple[str, ...] = BACKBONE_ATOMS
) -> dict[tuple, int]:
    mask = np.isin(array.chain_id, chains) & np.isin(array.atom_name, atoms)
    return {
        (
            str(array.chain_id[i]),
            int(array.res_id[i]),
            str(array.ins_code[i]),
            str(array.atom_name[i]),
        ): int(i)
        for i in np.where(mask)[0]
    }


def target_backbone_rmsd(
    array: AtomArray, reference: AtomArray, target_chains: list[str]
) -> float:
    """RMSD of the target backbone against the structure it was generated from.

    Computed **without superposition**, so a rigid-body displacement shows up
    here as a large number. That is informative but is *not* on its own a
    failure: RFdiffusion recentres its input, so its output is routinely
    translated by tens of angstroms from the structure it was given. Measured
    at 81 A RMS on the first NIM design of 9 Nov 2026.

    What must not happen is the target changing *shape*, and that is what
    :func:`superimpose_on_target` measures. Use this function to report the
    offset, that one to decide whether the geometry survived.

    Only atoms present in both structures are compared.

    Raises:
        ValueError: if the two share no backbone atoms at matching positions.
    """
    generated = _backbone_index(array, target_chains)
    original = _backbone_index(reference, target_chains)

    shared = sorted(set(generated) & set(original))
    if not shared:
        raise ValueError(
            "generated structure and reference share no target backbone atoms; "
            "chain ids or residue numbering differ"
        )

    a = array.coord[[generated[key] for key in shared]]
    b = reference.coord[[original[key] for key in shared]]
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def superimpose_on_target(
    array: AtomArray, reference: AtomArray, target_chains: list[str]
) -> tuple[AtomArray, float]:
    """Move a generated structure into the reference's frame, by its target.

    RFdiffusion holds the target rigid during denoising but writes the result
    recentred, so the output sits in a different frame from the stage-04 target
    it was conditioned on. Superimposing on the *target* backbone alone undoes
    exactly that rigid-body move and carries the binder along with it, which is
    legitimate precisely because the target was never allowed to flex.

    The returned RMSD is therefore a shape check: small when the target
    survived, and large when the output was misparsed, the contig was wrong, or
    the chains were reassigned. It is computed over :data:`FRAME_ATOMS` only,
    because the carbonyl oxygen is rebuilt rather than carried through.

    Returns:
        The whole array transformed into the reference frame, and the RMSD of
        the target backbone after superposition.

    Raises:
        ValueError: if the two share no backbone atoms at matching positions.
    """
    generated = _backbone_index(array, target_chains, FRAME_ATOMS)
    original = _backbone_index(reference, target_chains, FRAME_ATOMS)

    shared = sorted(set(generated) & set(original))
    if not shared:
        raise ValueError(
            "generated structure and reference share no target backbone atoms; "
            "chain ids or residue numbering differ"
        )

    mobile = array[[generated[key] for key in shared]]
    fixed = reference[[original[key] for key in shared]]

    fitted, transform = struc.superimpose(fixed, mobile)
    rmsd = float(struc.rmsd(fixed, fitted))
    return transform.apply(array), rmsd


def splice_binder(
    reference: AtomArray, generated: AtomArray, binder: str
) -> AtomArray:
    """Full-atom target plus generated binder, in the reference's frame.

    RFdiffusion emits backbone atoms only - for the target as well as the
    binder - so its output cannot be used directly for interface analysis: a
    side-chain-free target has nothing to bury, and Asp6's carboxylate, the
    whole point of this design problem, is simply absent.

    Since the target was held rigid, its side chains are known exactly: they
    are the ones in the stage-04 reference. This reattaches them by keeping the
    reference intact and adding only the generated chain. ``generated`` must
    already be in the reference frame - see :func:`superimpose_on_target`.

    The binder itself stays backbone-only, which is correct at this stage: it
    has no sequence yet.
    """
    chain = generated[generated.chain_id == binder]
    if chain.array_length() == 0:
        raise ValueError(f"generated structure has no chain {binder}")
    return reference + chain


def secondary_structure_fraction(array: AtomArray, chain_id: str) -> dict[str, float]:
    """Helix, strand and coil fractions of one chain.

    Uses biotite's P-SEA implementation, which needs only CA positions and so
    works on a sequence-free backbone. Reported because noise scale 0.5 is
    expected to bias strongly toward helical topologies, and that expectation
    should be checked rather than assumed.

    Raises:
        ValueError: if the chain has no residues.
    """
    chain = array[array.chain_id == chain_id]
    if chain.array_length() == 0:
        raise ValueError(f"chain {chain_id} has no atoms")

    sse = struc.annotate_sse(chain)
    if len(sse) == 0:
        raise ValueError(f"chain {chain_id} yielded no secondary-structure annotation")

    helix = float(np.count_nonzero(sse == SSE_HELIX)) / len(sse)
    strand = float(np.count_nonzero(sse == SSE_STRAND)) / len(sse)
    return {
        "helix": round(helix, 3),
        "strand": round(strand, 3),
        "coil": round(1.0 - helix - strand, 3),
    }


#: Empirical radius of gyration of a compact globular protein, Rg = A * N^B
#: angstroms for N residues. Used only to put Rg on a scale that is comparable
#: across binder lengths - the absolute constant matters far less than the ratio.
RG_PREFACTOR_A = 2.2
RG_EXPONENT = 0.38


def radius_of_gyration(array: AtomArray, chain_id: str) -> float:
    """Radius of gyration of one chain, over its CA atoms."""
    ca = array.coord[(array.chain_id == chain_id) & (array.atom_name == "CA")]
    if len(ca) == 0:
        raise ValueError(f"chain {chain_id} has no CA atoms")
    return float(np.sqrt(((ca - ca.mean(axis=0)) ** 2).sum(axis=1).mean()))


def compactness(array: AtomArray, chain_id: str) -> float:
    """Radius of gyration relative to a compact globule of the same length.

    1.0 is globular; a single extended alpha-helix comes out around 2.5, since
    70 residues of helix span roughly 100 angstroms.

    This exists because RFdiffusion at default noise produces many designs that
    are a single long helix. Those satisfy every interface criterion - a helix
    laid across the groove buries plenty of surface - while being unfoldable as
    standalone proteins, and ProteinMPNN answers them with low-complexity
    poly-alanine because there is no core to pack. The authors recommend
    noise_scale 0.5 for interfaces precisely because it improves designability,
    and the hosted API this project fell back to cannot set it.
    """
    ca_count = int(
        np.count_nonzero((array.chain_id == chain_id) & (array.atom_name == "CA"))
    )
    if ca_count == 0:
        raise ValueError(f"chain {chain_id} has no CA atoms")
    ideal = RG_PREFACTOR_A * ca_count**RG_EXPONENT
    return radius_of_gyration(array, chain_id) / ideal


def chain_length(array: AtomArray, chain_id: str) -> int:
    """Number of residues in one chain."""
    chain = array[array.chain_id == chain_id]
    if chain.array_length() == 0:
        raise ValueError(f"chain {chain_id} has no atoms")
    return struc.get_residue_count(chain)
