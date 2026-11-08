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


def _backbone_index(array: AtomArray, chains: list[str]) -> dict[tuple, int]:
    mask = np.isin(array.chain_id, chains) & np.isin(array.atom_name, BACKBONE_ATOMS)
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

    Computed **without superposition**, deliberately. RFdiffusion holds the
    target fixed at every denoising step, so the coordinates should be identical
    to the input. A superimposed RMSD would hide a rigid-body shift, which is
    exactly the failure this check exists to catch - it would mean the output
    was misparsed, the contig was wrong, or the chains were reassigned.

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


def chain_length(array: AtomArray, chain_id: str) -> int:
    """Number of residues in one chain."""
    chain = array[array.chain_id == chain_id]
    if chain.array_length() == 0:
        raise ValueError(f"chain {chain_id} has no atoms")
    return struc.get_residue_count(chain)
