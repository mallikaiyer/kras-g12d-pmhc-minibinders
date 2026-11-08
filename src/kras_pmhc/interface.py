"""Interface residue identification by relative solvent accessibility.

Implements the classification of Levy, *J Mol Biol* 403:660 (2010): a residue is
assigned to one of five classes from its relative SASA before and after complex
formation. This is the project-wide definition of an interface residue, and it
replaces distance cutoffs everywhere one was previously used.

Why rSASA rather than a distance cutoff: a cutoff counts residues that are near
but not occluded, and says nothing about burial. The classification separates
``core`` (buried on binding) from ``rim`` (at the edge, still solvated), which is
exactly the distinction this project's mechanistic claim rests on - that a
designed reader group is *buried and desolvated* against the aspartate, not
merely nearby.

The relative scale uses the Tien et al. (2013) theoretical maxima from
``constants.MAX_ASA_TIEN_2013``. Those were derived with their own radii and
probe, so mixing them with biotite's ProtOr radii introduces a small systematic
offset in rSASA. That is inherent to using tabulated maxima and is acceptable
here because every comparison in this project is made between structures treated
identically; the absolute rSASA value is never interpreted on its own.
"""

from __future__ import annotations

from dataclasses import dataclass

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray

from .constants import MAX_ASA_TIEN_2013

#: Relative SASA above which a residue counts as solvent-exposed (Levy 2010).
RSASA_THRESHOLD = 0.25

#: Absolute SASA loss, A^2, above which a residue counts as changed on binding.
#: A tolerance rather than ``> 0`` because Shrake-Rupley is a numerical
#: quadrature: identical geometry can differ by a tiny amount between the
#: monomer and complex evaluations, and that noise must not create interface
#: residues out of nothing.
DEFAULT_TOLERANCE_A2 = 0.1

#: Shrake-Rupley sampling, fixed so results are reproducible and reportable.
#: biotite's default is 100 points; 500 reduces the quadrature noise that
#: ``DEFAULT_TOLERANCE_A2`` exists to absorb.
SASA_POINT_NUMBER = 500
SASA_PROBE_RADIUS = 1.4
SASA_VDW_RADII = "ProtOr"

#: The three classes that together make up an interface (Levy 2010).
INTERFACE_CLASSES = frozenset({"support", "rim", "core"})

#: The two that are buried by the partner rather than remaining solvated.
BURIED_INTERFACE_CLASSES = frozenset({"support", "core"})


@dataclass(frozen=True)
class ResidueClass:
    """One residue's accessibility before and after complex formation.

    ``sasa_monomer`` is measured with the partner group's atoms removed at
    unchanged coordinates, so the difference reflects occlusion alone and not
    any conformational response.
    """

    chain_id: str
    res_id: int
    ins_code: str
    res_name: str
    group: str
    sasa_monomer: float
    sasa_complex: float
    rsasa_monomer: float
    rsasa_complex: float
    delta_rsasa: float
    changed: bool
    klass: str

    @property
    def sasa_loss(self) -> float:
        return self.sasa_monomer - self.sasa_complex


def per_residue_sasa(array: AtomArray) -> tuple[list[tuple], np.ndarray]:
    """Sum Shrake-Rupley atom areas within each residue.

    Returns:
        ``(keys, areas)`` where each key is
        ``(chain_id, res_id, ins_code, res_name)`` in the order the residues
        appear in ``array``, and ``areas`` is the matching area per residue.

    Raises:
        ValueError: if ``array`` holds no atoms.
    """
    if array.array_length() == 0:
        raise ValueError("cannot compute SASA of an empty atom array")

    per_atom = struc.sasa(
        array,
        probe_radius=SASA_PROBE_RADIUS,
        point_number=SASA_POINT_NUMBER,
        vdw_radii=SASA_VDW_RADII,
    )
    # nansum, not sum: ProtOr has no radius for some atoms and sasa returns NaN
    # for them. Treating those as zero area is the only sane choice, and it is
    # why the heavy-atom-only structures this project builds are preferred.
    areas = struc.apply_residue_wise(array, per_atom, np.nansum)

    starts = struc.get_residue_starts(array)
    keys = [
        (
            str(array.chain_id[i]),
            int(array.res_id[i]),
            str(array.ins_code[i]),
            str(array.res_name[i]),
        )
        for i in starts
    ]
    return keys, areas


def _validate_groups(
    array: AtomArray, group_a: list[str], group_b: list[str]
) -> None:
    if not group_a or not group_b:
        raise ValueError("both groups must name at least one chain")

    overlap = set(group_a) & set(group_b)
    if overlap:
        raise ValueError(f"groups overlap on chain(s) {sorted(overlap)}")

    present = set(np.unique(array.chain_id).tolist())
    missing = (set(group_a) | set(group_b)) - present
    if missing:
        raise ValueError(
            f"chain(s) {sorted(missing)} not in structure; present: {sorted(present)}"
        )


def _relative(sasa: float, res_name: str) -> float:
    try:
        return sasa / MAX_ASA_TIEN_2013[res_name]
    except KeyError as exc:
        raise KeyError(
            f"no maximum ASA tabulated for residue {res_name!r}; "
            "Tien et al. 2013 covers the 20 standard amino acids only"
        ) from exc


def _assign_class(
    changed: bool, rsasa_monomer: float, rsasa_complex: float, threshold: float
) -> str:
    """Levy's decision tree, in his order.

    ``support`` is tested before ``rim`` and ``core`` deliberately: a residue
    already buried within its own subunit is support whatever happens to its
    small remaining accessibility.
    """
    if not changed:
        return "interior" if rsasa_monomer < threshold else "surface"
    if rsasa_monomer < threshold:
        return "support"
    return "rim" if rsasa_complex >= threshold else "core"


def classify(
    array: AtomArray,
    group_a: list[str],
    group_b: list[str],
    tolerance: float = DEFAULT_TOLERANCE_A2,
    threshold: float = RSASA_THRESHOLD,
) -> list[ResidueClass]:
    """Classify every residue of two chain groups across their mutual interface.

    Two *groups* of chains rather than two chains: the literal per-chain reading
    of the method counts every internal interface, so on a binder plus
    HLA plus peptide system the HLA-peptide interface would swamp the signal
    being measured. Single-chain groups reproduce the per-chain case exactly.

    Only chains named in ``group_a`` or ``group_b`` take part; any other chain in
    ``array`` is excluded from both the complex and the monomer calculations, so
    the result depends on nothing that was not asked for.

    Args:
        array: structure holding at least the named chains.
        group_a: chain ids forming one side of the interface.
        group_b: chain ids forming the other side.
        tolerance: absolute SASA loss, A^2, above which a residue counts as
            changed on binding. Below it, ``delta_rsasa`` is reported as 0.
        threshold: relative SASA separating buried from exposed.

    Returns:
        One :class:`ResidueClass` per residue, ordered as in ``array``.

    Raises:
        ValueError: if a group is empty, the groups overlap, or a named chain is
            absent.
        KeyError: if a residue has no tabulated maximum ASA.
    """
    _validate_groups(array, group_a, group_b)

    in_a = np.isin(array.chain_id, group_a)
    in_b = np.isin(array.chain_id, group_b)
    complex_array = array[in_a | in_b]

    complex_keys, complex_areas = per_residue_sasa(complex_array)
    complex_sasa = dict(zip(complex_keys, complex_areas, strict=True))

    records: list[ResidueClass] = []
    for group_name, chains in (("a", group_a), ("b", group_b)):
        monomer = array[np.isin(array.chain_id, chains)]
        keys, areas = per_residue_sasa(monomer)

        for key, sasa_monomer in zip(keys, areas, strict=True):
            chain_id, res_id, ins_code, res_name = key
            sasa_complex = float(complex_sasa[key])
            sasa_monomer = float(sasa_monomer)

            rsasa_monomer = _relative(sasa_monomer, res_name)
            rsasa_complex = _relative(sasa_complex, res_name)

            changed = (sasa_monomer - sasa_complex) > tolerance
            delta_rsasa = rsasa_monomer - rsasa_complex if changed else 0.0

            records.append(
                ResidueClass(
                    chain_id=chain_id,
                    res_id=res_id,
                    ins_code=ins_code,
                    res_name=res_name,
                    group=group_name,
                    sasa_monomer=round(sasa_monomer, 3),
                    sasa_complex=round(sasa_complex, 3),
                    rsasa_monomer=round(rsasa_monomer, 4),
                    rsasa_complex=round(rsasa_complex, 4),
                    delta_rsasa=round(delta_rsasa, 4),
                    changed=changed,
                    klass=_assign_class(
                        changed, rsasa_monomer, rsasa_complex, threshold
                    ),
                )
            )
    return records


def interface_residues(records: list[ResidueClass]) -> list[ResidueClass]:
    """Support, rim and core - the interface as Levy defines it."""
    return [r for r in records if r.klass in INTERFACE_CLASSES]


def buried_interface_residues(records: list[ResidueClass]) -> list[ResidueClass]:
    """Support and core - the part of the interface the partner desolvates."""
    return [r for r in records if r.klass in BURIED_INTERFACE_CLASSES]


def interface_area_by_chain(records: list[ResidueClass]) -> dict[str, float]:
    """Total SASA lost, A^2, by interface residues of each chain.

    Used to ask which chain of a multi-chain target a binder actually engages -
    a binder that only contacts the HLA helices and never the peptide is
    useless here, and this is the number that detects it.
    """
    totals: dict[str, float] = {}
    for record in interface_residues(records):
        totals[record.chain_id] = totals.get(record.chain_id, 0.0) + record.sasa_loss
    return {chain: round(area, 3) for chain, area in totals.items()}


def find(
    records: list[ResidueClass], chain_id: str, res_id: int
) -> ResidueClass | None:
    """The record for one residue, or ``None`` if it is not present."""
    for record in records:
        if record.chain_id == chain_id and record.res_id == res_id:
            return record
    return None
