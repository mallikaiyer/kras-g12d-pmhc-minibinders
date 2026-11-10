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


# --- interface analysis (Levy 2010; Tien et al. 2013) ---------------------
# Theoretical maximum solvent-accessible area per residue, A^2, from Tien
# et al., PLoS ONE 8:e80635 (2013), used as the denominator for relative SASA.
# These are the *theoretical* (Gly-X-Gly tripeptide) maxima rather than the
# empirical ones: the empirical values are drawn from observed structures and
# so already reflect burial, which is the quantity being measured here.
MAX_ASA_TIEN_2013 = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}


# --- challenge target sequences -------------------------------------------
# Taken from the Challenge 5 page, supplied by the user on 9 Nov 2026. These are
# the sequences the submission must name as the target, independent of the
# structures this project was built from.
#
# Verified rather than trusted: the first 180 residues of the heavy chain match
# chain A of data/targets_relaxed/7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb
# exactly, 180 of 180. That confirms both the allele and - more usefully - that
# the residue numbering used throughout stages 00-05 is the numbering the
# challenge uses. tests/test_constants.py pins the agreement so a later edit to
# either the sequence or the target cannot silently break it.
HLA_A1101_HEAVY_CHAIN = (
    "GSHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGT"
    "LRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHAAEQQRAYL"
    "EGRCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGT"
    "FQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWELSSQPT"
)

#: Residues of the heavy chain covered by the design target crop (alpha1/alpha2).
HLA_DESIGN_TARGET_RESIDUES = 180

# beta-2 microglobulin, from the same challenge page, supplied 9 Nov 2026.
# Verified against 7OW3 chain B: residues 1-99 of the structure match this
# exactly, 99 of 99. The structure additionally carries an initiator Met
# numbered 0, so the challenge specifies the *mature* chain and the crystal the
# expressed one - a distinction the stage-11 construct description has to get
# right. beta-2m is not part of the design target: stage 03 cropped it away
# along with alpha3, because the binder never sees it.
B2M_SEQUENCE = (
    "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYA"
    "CRVNHVTLSQPKIVKWDRDM"
)

#: The initiator Met present in the crystal structures but not in the mature
#: sequence above, numbered residue 0 in 7OW3 chain B.
B2M_STRUCTURE_HAS_INITIATOR_MET = True


# --- sequence liability screens -------------------------------------------
# Kyte & Doolittle (1982) J Mol Biol 157:105. The competition's design-selection
# prompt names mean Kyte-Doolittle hydropathy as one of three screens applied to
# every design, alongside pLDDT and ipSAE, so it is reported per sequence.
# Positive is hydrophobic; a de novo binder averaging strongly positive is a
# solubility and aggregation risk.
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
