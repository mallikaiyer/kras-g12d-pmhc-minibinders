#!/usr/bin/env python
"""Stage 08 - score the cofolded complexes.

Turns each predicted structure into the three numbers the project's ranking
needs, in its stated order of precedence:

1. **Mutation specificity.** The same binder is folded against the G12D and
   wild-type peptides, which differ at position 6 alone. Specificity is the
   difference between the two - in PRODIGY binding free energy, in the Levy
   burial of residue 6, and in interface confidence. A binder that does equally
   well on both reads the groove, not the neoepitope.
2. **Peptide-dependent binding.** The peptide's share of the target-side
   interface in the G12D complex, by the Levy (2010) rSASA classification used
   everywhere in this project.
3. **Affinity.** PRODIGY's predicted dG for the G12D complex. Last, and only
   meaningful for designs that clear the first two.

This is the stage that cofolding made possible: PRODIGY counts atomic contacts
across an interface, so it needs side chains on both partners, and stage 05
produced backbones only.

Interfaces are identified by rSASA, never by a distance cutoff. biotite does all
of this project's own parsing; PRODIGY requires a Bio.PDB model, so a converted
copy is handed to it at that boundary and nowhere else.

Usage:
    python scripts/08_score_cofolded.py
    python scripts/08_score_cofolded.py --results data/cofolding_nim/RESULTS.tsv
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import warnings
from pathlib import Path

import biotite.structure as struc
from biotite.sequence import ProteinSequence
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import backbones, interface  # noqa: E402
from kras_pmhc.constants import MUTANT_POSITION  # noqa: E402

DEFAULT_RESULTS = ROOT / "data" / "cofolding_nim" / "RESULTS.tsv"
OUT = ROOT / "results" / "cofolding" / "SCORES.tsv"

# Chain ids in the *output* are not the ids that were submitted. Boltz-2 was
# sent polymers with ids A (HLA), C (peptide), B (binder) and returned them
# relabelled sequentially in submission order: A, B, C. Trusting the submitted
# ids put the binder where the peptide was expected, and the only reason that
# did not become a table of binder-against-itself interface areas is that
# EXPECTED_RESIDUE_6 below checks the identity of residue 6.
#
# So chains are resolved by sequence, not by label. Lengths alone would be
# ambiguous if a binder were ever 10 or 180 residues long.
HLA_LENGTH = 180
PEPTIDE_LENGTH = 10

#: Residue expected at peptide position 6 in each condition. Checked rather than
#: assumed: if these are the wrong way round, every specificity number inverts.
EXPECTED_RESIDUE_6 = {"g12d": "ASP", "wt": "GLY"}

#: pair_chains_iptm is keyed by *submission* index, which is unaffected by the
#: relabelling: 0 HLA, 1 peptide, 2 binder.
IPTM_HLA, IPTM_PEPTIDE, IPTM_BINDER = "0", "1", "2"

FIELDS = (
    "job_id", "backbone_id", "sample", "condition",
    "residue6", "residue6_class", "residue6_sasa_loss", "residue6_delta_rsasa",
    "peptide_interface_area", "hla_interface_area", "peptide_share",
    "n_interface_binder", "buried_interface_area",
    "binder_compactness", "binder_helix",
    "iptm_binder_peptide", "iptm_binder_hla", "iptm_hla_peptide",
    "complex_plddt", "prodigy_dg", "prodigy_kd", "n_contacts",
)


def load_structure(path: Path):
    """Read a prediction, mmCIF or PDB, first model only.

    Hydrogens are dropped. Everything downstream - the Levy rSASA
    classification, PRODIGY's contact counting - is parameterised on heavy
    atoms, and the minimised structures from stage 10 carry hydrogens that
    would change SASA values and make minimised and unminimised scores
    incomparable.
    """
    if path.suffix == ".cif":
        array = get_structure(CIFFile.read(str(path)), model=1)
    else:
        array = PDBFile.read(str(path)).get_structure(model=1)
    return array[array.element != "H"]


def chain_sequence(array, chain_id: str) -> str:
    chain = array[array.chain_id == chain_id]
    starts = struc.get_residue_starts(chain)
    return "".join(
        ProteinSequence.convert_letter_3to1(name) for name in chain.res_name[starts]
    )


def resolve_chains(array, binder_sequence: str, peptide_sequence: str) -> dict:
    """Work out which output chain is which, from the sequences themselves.

    Raises:
        SystemExit: if any of the three cannot be identified unambiguously.
            Guessing here would silently compute the wrong interface.
    """
    found: dict[str, str] = {}
    for chain_id in sorted(set(array.chain_id)):
        sequence = chain_sequence(array, chain_id)
        if sequence == binder_sequence:
            found["binder"] = chain_id
        elif sequence == peptide_sequence:
            found["peptide"] = chain_id
        elif len(sequence) == HLA_LENGTH:
            found["hla"] = chain_id

    missing = {"hla", "peptide", "binder"} - set(found)
    if missing:
        present = {
            c: f"{len(chain_sequence(array, c))} aa"
            for c in sorted(set(array.chain_id))
        }
        raise SystemExit(
            f"could not identify {sorted(missing)} among chains {present}; "
            "the model does not contain the sequences that were submitted"
        )
    return found


def to_pdb_text(array) -> str:
    """mmCIF to PDB text, for handing to PRODIGY."""
    out = PDBFile()
    out.set_structure(array)
    handle = io.StringIO()
    out.write(handle)
    return handle.getvalue()


def prodigy_scores(pdb_text: str, name: str, selection: list[str]) -> dict:
    """PRODIGY's dG and contact count for binder against target.

    Imported here rather than at module scope so this script still runs, with
    these columns blank, in an environment without PRODIGY.
    """
    try:
        from Bio.PDB import PDBParser
        from prodigy_prot.modules.prodigy import Prodigy
    except ImportError:
        return {}

    parser = PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure(name, io.StringIO(pdb_text))
        model = structure[0]
        # Binder against both target chains, matching the Levy grouping used
        # for every interface number in this project. The selection is passed in
        # because the output chain ids are resolved per model, not assumed.
        prodigy = Prodigy(model, name=name, selection=selection)
        prodigy.predict()
    result = prodigy.as_dict()
    return {
        "prodigy_dg": round(result.get("ba_val", float("nan")), 3),
        "prodigy_kd": f"{result.get('kd', float('nan')):.3e}",
        "n_contacts": result.get("ICs_total", result.get("ic_total", "")),
    }


def score(path: Path, job: dict) -> dict:
    array = load_structure(path)
    condition = job["condition"]

    chains = resolve_chains(array, job["binder_sequence"], job["peptide_sequence"])
    binder_chain = chains["binder"]
    peptide_chain = chains["peptide"]
    hla_chain = chains["hla"]
    target_chains = [hla_chain, peptide_chain]

    records = interface.classify(array, [binder_chain], target_chains)
    residue6 = interface.find(records, peptide_chain, MUTANT_POSITION)
    if residue6 is None:
        raise SystemExit(f"{job['job_id']}: peptide residue {MUTANT_POSITION} absent")

    expected = EXPECTED_RESIDUE_6[condition]
    if residue6.res_name != expected:
        raise SystemExit(
            f"{job['job_id']}: peptide residue {MUTANT_POSITION} is "
            f"{residue6.res_name}, expected {expected} for condition "
            f"'{condition}' - the conditions are mislabelled and every "
            "specificity number would be inverted"
        )

    areas = interface.interface_area_by_chain(records)
    peptide = areas.get(peptide_chain, 0.0)
    hla = areas.get(hla_chain, 0.0)
    total = peptide + hla

    binder_records = [r for r in records if r.chain_id == binder_chain]
    sse = backbones.secondary_structure_fraction(array, binder_chain)

    pairs = json.loads(job["pair_chains_iptm"]) if job.get("pair_chains_iptm") else {}

    def pair_iptm(i: str, j: str):
        try:
            return round((pairs[i][j] + pairs[j][i]) / 2, 4)
        except (KeyError, TypeError):
            return ""

    row = {
        "job_id": job["job_id"],
        "backbone_id": job["backbone_id"],
        "sample": job["sample"],
        "condition": condition,
        "residue6": residue6.res_name,
        "residue6_class": residue6.klass,
        "residue6_sasa_loss": round(residue6.sasa_loss, 2),
        "residue6_delta_rsasa": residue6.delta_rsasa,
        "peptide_interface_area": round(peptide, 2),
        "hla_interface_area": round(hla, 2),
        "peptide_share": round(peptide / total, 4) if total else 0.0,
        "n_interface_binder": sum(
            1 for r in binder_records if r.klass in interface.INTERFACE_CLASSES
        ),
        "buried_interface_area": round(
            sum(r.sasa_loss for r in records if r.klass in interface.BURIED_INTERFACE_CLASSES),
            2,
        ),
        "binder_compactness": round(backbones.compactness(array, binder_chain), 3),
        "binder_helix": sse["helix"],
        "iptm_binder_peptide": pair_iptm(IPTM_BINDER, IPTM_PEPTIDE),
        "iptm_binder_hla": pair_iptm(IPTM_BINDER, IPTM_HLA),
        "iptm_hla_peptide": pair_iptm(IPTM_HLA, IPTM_PEPTIDE),
        "complex_plddt": job.get("complex_plddt", ""),
        "prodigy_dg": "",
        "prodigy_kd": "",
        "n_contacts": "",
    }
    row.update(prodigy_scores(
        to_pdb_text(array),
        job["job_id"],
        [binder_chain, f"{hla_chain},{peptide_chain}"],
    ))
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument(
        "--structure",
        default="model.cif",
        help="file to score inside each job directory: model.cif (the raw "
             "prediction, default) or model_min.pdb (minimised by stage 10)",
    )
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    if not args.results.is_file():
        raise SystemExit(f"{args.results} not found - run scripts/07b_nim_cofolding.py")
    with args.results.open() as handle:
        jobs = list(csv.DictReader(handle, delimiter="\t"))

    root = args.results.parent
    rows = []
    print(f"scoring {len(jobs)} cofolded complex(es) from {root}")
    print(f"structure: {args.structure}\n")
    for job in jobs:
        model = root / job["job_id"] / args.structure
        if not model.is_file():
            print(f"  {job['job_id']:<22} missing {args.structure}")
            continue
        row = score(model, job)
        rows.append(row)
        print(f"  {row['job_id']:<22} res6 {row['residue6']} {row['residue6_class']:<8}"
              f" loss {float(row['residue6_sasa_loss']):>6.1f}"
              f" share {float(row['peptide_share']):>5.0%}"
              f" dG {row['prodigy_dg'] or 'n/a':>7}")

    if not rows:
        raise SystemExit("nothing scored")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nwrote {args.out.relative_to(ROOT)}  (from {args.structure})")
    if not any(row["prodigy_dg"] != "" for row in rows):
        print("PRODIGY columns are blank - install prodigy-prot to fill them")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
