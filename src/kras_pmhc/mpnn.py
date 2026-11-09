"""ProteinMPNN sequence design via the hosted NVIDIA NIM endpoint.

Stage 06. RFdiffusion emits a poly-glycine backbone, so a sequence has to be
designed onto it before anything can be expressed, folded or scored. The binder
chain is designed while the HLA and the peptide are held as fixed context -
ProteinMPNN is conditioned on the target's sequence as well as its coordinates,
so it knows peptide position 6 is an aspartate even though it cannot see which
way the carboxylate points.

Two measured facts shape the parameters, both established by probing the
endpoint before any sequence was kept:

* ``use_soluble_model`` is *worse* here. On nim3_6 it returned sequences that
  were 71-80% hydrophobic with runs of 43-53 alanines, against 43-63% for the
  default model. It is trained on soluble monomers and this is an interface.
* ``random_seed`` is **not honoured**. Two identical requests with seed 0
  returned different sequences, so stage 06 is not reproducible by re-running
  it. The raw multi-FASTA response is therefore preserved verbatim: it is the
  only record of what was designed, and regenerating would silently produce
  different sequences.

The first record of a ProteinMPNN response is always the *input* sequence - for
an RFdiffusion backbone, a run of glycines - and must be skipped rather than
counted as a design.
"""

from __future__ import annotations

import re
import urllib.request
from dataclasses import dataclass

from kras_pmhc import nim

#: Sampling temperature. 0.1 is the conservative end of the 0.1-0.3 the NIM
#: documentation recommends; higher values trade confidence for diversity, and
#: there are only nine backbones to work with.
DEFAULT_TEMPERATURE = 0.1

#: Sequences per backbone.
DEFAULT_NUM_SEQUENCES = 10

#: Header of the record carrying the input sequence rather than a design.
INPUT_RECORD_PREFIX = ">input"

#: Fields ProteinMPNN reports per designed sequence.
_HEADER_FIELD = re.compile(r"(\w+)=([-\d.eE]+)")


class MpnnError(RuntimeError):
    """A response that cannot be interpreted as sequence designs."""


@dataclass(frozen=True)
class Sequence:
    """One designed sequence and the scores ProteinMPNN reported for it."""

    sample: int
    temperature: float
    sequence: str
    score: float
    global_score: float
    seq_recovery: float

    @property
    def length(self) -> int:
        return len(self.sequence)


def build_payload(
    pdb_text: str,
    design_chains: list[str],
    num_sequences: int = DEFAULT_NUM_SEQUENCES,
    temperature: float = DEFAULT_TEMPERATURE,
    seed: int | None = None,
) -> dict:
    """Assemble one request.

    ``input_pdb_chains`` names the chains to *design*; everything else in the
    structure is kept as fixed context, which the response confirms by echoing
    ``fixed_chains`` and ``designed_chains``.
    """
    if not design_chains:
        raise MpnnError("no chain to design")
    if not 1 <= num_sequences <= 100:
        raise MpnnError(f"num_sequences={num_sequences} outside the documented 1-100")
    if not 0 < temperature <= 1:
        raise MpnnError(f"temperature={temperature} outside the documented (0, 1]")

    payload = {
        "input_pdb": nim.atom_lines(pdb_text),
        "input_pdb_chains": list(design_chains),
        "num_seq_per_target": num_sequences,
        "sampling_temp": [temperature],
        "ca_only": False,
        # Measured, not assumed: the soluble model gave 71-80% hydrophobic
        # sequences here against 43-63% for the default.
        "use_soluble_model": False,
    }
    if seed is not None:
        # Sent for the record even though the service ignores it - see the
        # module docstring. Its presence in the log is a claim about what was
        # requested, not about reproducibility.
        payload["random_seed"] = seed
    return payload


def parse_mfasta(mfasta: str, expected_length: int | None = None) -> list[Sequence]:
    """Pull the designed sequences out of a ProteinMPNN multi-FASTA response.

    Skips the leading ``>input`` record, which holds the backbone's own
    sequence - all glycine, for RFdiffusion output.

    Raises:
        MpnnError: if no designs are present, or one is not the expected
            length. A length mismatch would mean the wrong chain was designed,
            and silently accepting it would attach a sequence to the wrong
            residues.
    """
    lines = [line.strip() for line in mfasta.splitlines() if line.strip()]
    designs: list[Sequence] = []

    for index, line in enumerate(lines):
        if not line.startswith(">") or line.startswith(INPUT_RECORD_PREFIX):
            continue
        if index + 1 >= len(lines) or lines[index + 1].startswith(">"):
            raise MpnnError(f"header with no sequence after it: {line[:80]}")

        fields = dict(_HEADER_FIELD.findall(line))
        sequence = lines[index + 1]
        if expected_length is not None and len(sequence) != expected_length:
            raise MpnnError(
                f"designed {len(sequence)} residues, expected {expected_length}; "
                "the wrong chain was designed"
            )
        designs.append(
            Sequence(
                sample=int(float(fields.get("sample", len(designs) + 1))),
                temperature=float(fields.get("T", DEFAULT_TEMPERATURE)),
                sequence=sequence,
                score=float(fields.get("score", "nan")),
                global_score=float(fields.get("global_score", "nan")),
                seq_recovery=float(fields.get("seq_recovery", "nan")),
            )
        )

    if not designs:
        raise MpnnError("response contained no designed sequences")
    return designs


def design(
    key: str,
    pdb_text: str,
    design_chains: list[str],
    num_sequences: int = DEFAULT_NUM_SEQUENCES,
    temperature: float = DEFAULT_TEMPERATURE,
    seed: int | None = None,
    expected_length: int | None = None,
    opener=urllib.request.urlopen,
) -> tuple[list[Sequence], str]:
    """Design sequences for one backbone.

    Returns the parsed designs and the raw multi-FASTA, which the caller is
    expected to keep: the service ignores the seed, so this response cannot be
    reproduced.
    """
    payload = build_payload(
        pdb_text, design_chains, num_sequences, temperature, seed
    )
    body = nim.post(nim.PROTEINMPNN_URL, key, payload, opener)

    mfasta = body.get("mfasta")
    if not mfasta:
        raise MpnnError(f"200 OK but no mfasta in response: {list(body)}")
    return parse_mfasta(mfasta, expected_length), mfasta
