"""Loading and chain-level interpretation of the pMHC structures.

The structures in this project have between one and four copies of the complex
in the asymmetric unit, and the chain naming differs between entries. Nothing
here hard-codes a chain letter: roles are assigned from residue counts and
from the presence of the HLA groove residues, so a chain mix-up fails loudly
rather than producing a plausible wrong number.
"""

from __future__ import annotations

from pathlib import Path

import biotite.structure as struc
import numpy as np
from biotite.structure import AtomArray
from biotite.structure.io.pdbx import CIFFile, get_structure

from .constants import HLA_ARG114, HLA_GLN70, STRUCTURES, THREE_LETTER

#: Default location of fetched structures, relative to the repository root.
DEFAULT_PDB_DIR = Path(__file__).resolve().parents[2] / "data" / "external" / "pdb"

# Chain-role thresholds, in residues. Class I heavy chains run ~275-300
# modelled residues, beta-2-microglobulin ~99, and the peptide is a 10-mer.
# TCR variable+constant domains overlap the heavy-chain range, which is why
# heavy-chain identification also requires the groove residues (see below).
MAX_PEPTIDE_RESIDUES = 15
B2M_RESIDUE_RANGE = (85, 125)
MIN_HEAVY_CHAIN_RESIDUES = 250

#: A peptide and its presenting groove are in van der Waals contact.
PEPTIDE_CONTACT_CUTOFF = 5.0


def load(pdb_id: str, pdb_dir: Path | None = None) -> AtomArray:
    """Read one structure as an ``AtomArray`` of protein atoms, model 1.

    Waters, ions and ligands are dropped. ``b_factor`` and ``occupancy`` are
    requested explicitly because biotite does not load them by default, and
    both are needed to judge how well-determined a residue actually is.
    Alternate conformations are reduced to the first, so downstream geometry
    sees exactly one coordinate per atom.
    """
    pdb_dir = Path(pdb_dir) if pdb_dir is not None else DEFAULT_PDB_DIR
    path = pdb_dir / f"{pdb_id.lower()}.cif"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `make fetch` (scripts/00_fetch_structures.py) first."
        )

    array = get_structure(
        CIFFile.read(str(path)),
        model=1,
        extra_fields=["b_factor", "occupancy"],
        altloc="first",
    )
    return array[struc.filter_amino_acids(array) & ~array.hetero]


def chain_residue_counts(array: AtomArray) -> dict[str, int]:
    """Number of distinct modelled residues per chain."""
    counts: dict[str, int] = {}
    for chain_id in np.unique(array.chain_id):
        chain = array[array.chain_id == chain_id]
        counts[str(chain_id)] = len(np.unique(chain.res_id))
    return counts


def _has_groove_residues(array: AtomArray, chain_id: str) -> bool:
    """Whether a chain carries the HLA groove markers Arg114 and Gln70.

    Used to separate the class I heavy chain from a TCR beta chain of similar
    length. Both markers must match by residue name, not merely by number.
    """
    chain = array[array.chain_id == chain_id]
    for res_id, expected in ((HLA_ARG114, "ARG"), (HLA_GLN70, "GLN")):
        names = np.unique(chain.res_name[chain.res_id == res_id])
        if expected not in names:
            return False
    return True


def peptide_chains(array: AtomArray) -> list[str]:
    """Chain IDs of the presented peptides, in file order."""
    counts = chain_residue_counts(array)
    return [c for c, n in counts.items() if n <= MAX_PEPTIDE_RESIDUES]


def heavy_chains(array: AtomArray) -> list[str]:
    """Chain IDs of class I heavy chains, identified by length and groove residues."""
    counts = chain_residue_counts(array)
    return [
        c
        for c, n in counts.items()
        if n >= MIN_HEAVY_CHAIN_RESIDUES and _has_groove_residues(array, c)
    ]


def pair_peptides_with_heavy_chains(array: AtomArray) -> dict[str, str]:
    """Map each peptide chain to the heavy chain presenting it.

    Pairing is by closest approach rather than by chain-letter convention,
    because the asymmetric unit ordering is not guaranteed.

    Raises:
        ValueError: if a peptide has no heavy chain in contact, which means
            either a parsing problem or a structure this code does not
            understand. Guessing would be worse than failing.
    """
    pairs: dict[str, str] = {}
    candidates = heavy_chains(array)
    if not candidates:
        raise ValueError("no class I heavy chain found (expected Arg114 and Gln70)")

    for pep in peptide_chains(array):
        pep_coords = array.coord[array.chain_id == pep]
        best, best_dist = None, np.inf
        for heavy in candidates:
            heavy_coords = array.coord[array.chain_id == heavy]
            d = float(
                np.linalg.norm(
                    pep_coords[:, None, :] - heavy_coords[None, :, :], axis=-1
                ).min()
            )
            if d < best_dist:
                best, best_dist = heavy, d
        if best is None or best_dist > PEPTIDE_CONTACT_CUTOFF:
            raise ValueError(
                f"peptide chain {pep} has no heavy chain within "
                f"{PEPTIDE_CONTACT_CUTOFF} A (closest {best_dist:.2f} A)"
            )
        pairs[pep] = best
    return pairs


def residue(array: AtomArray, chain_id: str, res_id: int) -> AtomArray:
    """All atoms of one residue. Empty if the residue is not modelled."""
    return array[(array.chain_id == chain_id) & (array.res_id == res_id)]


def atom_coord(
    array: AtomArray, chain_id: str, res_id: int, atom_name: str
) -> np.ndarray | None:
    """Coordinates of one named atom, or ``None`` if it is not modelled.

    ``None`` means "not observed in the crystal", which is a routine and
    informative outcome here - several peptide copies lack the central
    residues entirely - so callers must handle it rather than treat it as an
    error.
    """
    sel = residue(array, chain_id, res_id)
    hit = sel.coord[sel.atom_name == atom_name]
    if len(hit) == 0:
        return None
    return np.asarray(hit[0], dtype=float)


def observed_sequence(array: AtomArray, chain_id: str, length: int) -> str:
    """One-letter sequence of a peptide chain, with ``-`` for unmodelled positions.

    Positions are read from author numbering 1..length, which is how the
    peptide is numbered in all five entries.
    """
    reverse = {v: k for k, v in THREE_LETTER.items()}
    out = []
    for res_id in range(1, length + 1):
        sel = residue(array, chain_id, res_id)
        if len(sel) == 0:
            out.append("-")
        else:
            out.append(reverse.get(str(sel.res_name[0]), "X"))
    return "".join(out)


def check_peptide_identity(array: AtomArray, pdb_id: str) -> dict[str, str]:
    """Verify every peptide copy matches the expected sequence where modelled.

    Returns a mapping of chain ID to observed sequence.

    Raises:
        ValueError: on any mismatch at a modelled position. This is the guard
            against analysing the wrong chain or the wrong entry.
    """
    expected = STRUCTURES[pdb_id.upper()].peptide
    observed: dict[str, str] = {}
    for pep in peptide_chains(array):
        seq = observed_sequence(array, pep, len(expected))
        observed[pep] = seq
        for position, (got, want) in enumerate(zip(seq, expected, strict=True), start=1):
            if got != "-" and got != want:
                raise ValueError(
                    f"{pdb_id} chain {pep}: position {position} is {got}, expected {want} "
                    f"(observed {seq}, expected {expected})"
                )
    if not observed:
        raise ValueError(f"{pdb_id}: no peptide chain found")
    return observed
