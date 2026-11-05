"""Target definitions and the literature values this project is calibrated against.

Everything here is a fact about the system, taken from the challenge page or
from Poole et al., Nat Commun 13:5333 (2022), doi:10.1038/s41467-022-32811-1.
Keeping them in one place means the analysis scripts assert against named
constants rather than against magic numbers scattered through the code.
"""

from __future__ import annotations

# --- the peptides ---------------------------------------------------------
# KRAS residues 7-16 (UniProt P01116). The two peptides differ only at
# position 6: Gly12 -> Asp12 in KRAS numbering.
PEPTIDE_G12D = "VVVGADGVGK"
PEPTIDE_WT = "VVVGAGGVGK"

#: The single position that distinguishes the two assay targets.
MUTANT_POSITION = 6

# Three-letter codes, for asserting against structure files.
THREE_LETTER = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}


# --- the structures -------------------------------------------------------
class Structure:
    """One PDB entry relevant to this challenge."""

    def __init__(self, pdb_id: str, peptide: str, has_tcr: bool, description: str):
        self.pdb_id = pdb_id
        self.peptide = peptide
        self.has_tcr = has_tcr
        self.description = description

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Structure({self.pdb_id})"


#: All five structures from Poole et al. 2022. Note that 7OW3 and 7PB2 are not
#: mentioned on the challenge page but are directly relevant: 7OW3 is the only
#: entry with modelled coordinates at peptide position 6 in a *free* complex.
STRUCTURES = {
    "7OW3": Structure("7OW3", PEPTIDE_WT, False, "HLA-A*11:01 + KRAS-WT, free (4 copies)"),
    "7OW4": Structure("7OW4", PEPTIDE_G12D, False, "HLA-A*11:01 + KRAS-G12D, free (4 copies)"),
    "7OW5": Structure("7OW5", PEPTIDE_WT, True, "HLA-A*11:01 + KRAS-WT + JDIa41b1 TCR"),
    "7OW6": Structure("7OW6", PEPTIDE_G12D, True, "HLA-A*11:01 + KRAS-G12D + JDIa41b1 TCR"),
    "7PB2": Structure("7PB2", PEPTIDE_G12D, True, "HLA-A*11:01 + KRAS-G12D + parent JDI TCR"),
}


# --- HLA residues of interest --------------------------------------------
# PDB numbering, as used by the challenge page and the paper. Arg114 forms the
# salt bridge to Asp6 when the peptide is in the closed (TCR-bound) state;
# Asp116 stabilises that arginine; Gln70 donates a hydrogen bond.
HLA_ARG114 = 114
HLA_GLN70 = 70
HLA_ASP116 = 116

#: HLA positions whose distance to peptide position 6 the survey reports.
#: Ala158 and Lys146 are included as negative controls: the paper places them
#: near the TCR interface but far from the mutation.
HLA_SURVEY_RESIDUES = {
    70: "GLN",
    114: "ARG",
    116: "ASP",
    146: "LYS",
    155: "GLN",
    158: "ALA",
}

#: Reference atom for "is position 6 pointing into the groove?". CZ is the
#: central carbon of the Arg114 guanidinium, i.e. the tip of the side chain
#: that Asp6 pairs with in the closed state.
GROOVE_REFERENCE_ATOM = ("ARG", HLA_ARG114, "CZ")


# --- TCR residues carrying the discrimination -----------------------------
# Per-residue MMPBSA decomposition from Poole et al., Fig. 3b-d and
# Supplementary Table S3. Chain D is TCR alpha and chain E is TCR beta in
# 7OW6; the survey asserts the residue identities to catch a chain mix-up.
#
# The point of this table is that the two largest contributors sit >13 A from
# the Asp6 carboxylate: the discrimination is long-range electrostatic, not
# contact-mediated. The paper's phrase "directly positioned above the mutant
# residue" refers to the CDR3-beta loop as a whole, not to these residues
# individually, and K70a is explicitly described as further away.
#
# Keyed by PDB entry, and deliberately populated for 7OW6 only: the
# decomposition was performed on the JDIa41b1 complex. 7PB2 contains the
# *parent* JDI TCR, a different protein whose chain assignment and residue
# numbering do not correspond - applying this table to it silently produces
# distances for unrelated residues.
TCR_DISCRIMINATING_RESIDUES = {
    # pdb_id: [(chain, res_id, expected_resname, label, ddG_kcal_per_mol or None)]
    "7OW6": [
        ("D", 50, "ARG", "R50a", 2.0),
        ("D", 70, "LYS", "K70a", 0.88),
        ("E", 94, "LYS", "K94b", None),
        ("E", 99, "GLN", "Q99b", None),
        ("E", 100, "HIS", "H100b", None),
        ("E", 101, "ASN", "N101b", None),
    ],
}

#: Atoms making up the Asp6 carboxylate, the group that would be read directly
#: by a designed binder. CB is excluded: it carries no charge and including it
#: flatters the apparent proximity of distant residues.
ASP_CARBOXYLATE_ATOMS = ("CG", "OD1", "OD2")


# --- experimental calibration (Poole et al. Tables 1-2) -------------------
# Affinity of the parent and matured TCRs. The key observation for design is
# that the parent binds G12D only weakly yet shows no detectable binding to
# wild-type, while affinity maturation progressively *destroyed* that
# selectivity. Micromolar affinity is therefore the target, not nanomolar.
TCR_AFFINITIES = [
    # (name, KD_G12D_molar, KD_WT_molar or None for no binding, window)
    ("JDI (parent)", 6.30e-05, None, None),
    ("JDIa9bwt (1st gen)", 2.61e-07, 2.58e-05, 98),
    ("JDIa41b1 (2nd gen)", 7.43e-10, 3.00e-06, 4000),
    ("JDIa96b35 (3rd gen)", 6.37e-11, 4.25e-07, 6000),
]

#: Isothermal titration / van't Hoff thermodynamics for JDIa41b1, kcal/mol.
#: The G12D complex is far more enthalpically favourable and pays a much
#: larger entropic penalty - the signature of an ordering transition on
#: binding, which is what the conformational-selection argument assumes.
THERMODYNAMICS = {
    "WT": {"dH": -9.56, "TdS": -1.82, "dG": -7.65, "dCp": -0.31},
    "G12D": {"dH": -20.32, "TdS": -8.13, "dG": -12.43, "dCp": -0.41},
}
