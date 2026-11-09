#!/usr/bin/env python
"""Stage 05 (fallback) - generate binder backbones via the NVIDIA NIM API.

Replaces the Colab notebook when no GPU is available. The output layout is
deliberately identical to the notebook's, so scripts/05_ingest_backbones.py
scores these backbones with no change at all:

    <out-dir>/backbones/*.pdb
    <out-dir>/run_log.tsv

The recorded cost of this route is that the API has no noise-scale field, so
designs are sampled at RFdiffusion's default noise rather than the 0.5 the
authors recommend for protein-protein interfaces. ``noise_scale`` is logged as
``default(api)`` rather than left blank, so no later stage can mistake these
backbones for notebook ones.

Usage:
    # smoke test: one design, prove the contigs and target are accepted
    python scripts/05b_nim_generate.py --n 1

    # pilot
    python scripts/05b_nim_generate.py --n 8

    # production to a wall-clock cutoff
    python scripts/05b_nim_generate.py --until 2026-11-09T22:00:00
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from datetime import UTC, datetime
from pathlib import Path

from biotite.structure.io.pdb import PDBFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import nim  # noqa: E402

DEFAULT_TARGET = "7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb"
RELAXED = ROOT / "data" / "targets_relaxed"
DEFAULT_OUT = ROOT / "data" / "nim_stage05"
DEFAULT_KEY_FILE = ROOT / ".nim_api_key"

#: Same values as notebooks/05_rfdiffusion_binder.ipynb cell 1, so the two
#: routes differ in transport and noise scale only.
CONTIGS = "A1-180/0 C1-10/0 60-75"
HOTSPOTS = ["C5", "C6", "C7"]
DIFFUSION_STEPS = 50

#: Matches the notebook's LOG_FIELDS exactly; the ingest script keys on
#: ``filename`` and reads ``batch`` and ``seed``.
LOG_FIELDS = (
    "batch", "seed", "design", "filename", "contigs", "hotspots",
    "noise_scale", "diffuser_T", "checkpoint", "elapsed_s", "written_utc",
)

#: Recorded in place of a number because the API cannot set it. A blank would
#: read as "not recorded"; this reads as "not settable".
NOISE_SCALE = "default(api)"
CHECKPOINT = "nim-hosted(unspecified)"


def read_log(out_dir: Path) -> list[dict]:
    path = out_dir / "run_log.tsv"
    if not path.is_file():
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def append_log(out_dir: Path, row: dict) -> None:
    """Append after every design, not at the end of the run.

    A rate limit or a dropped connection halfway through must not cost the
    provenance of the backbones already written.
    """
    path = out_dir / "run_log.tsv"
    exists = path.is_file()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOG_FIELDS, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def next_seed(log: list[dict], seed_start: int) -> int:
    """Resume rather than repeat: a re-run continues the seed sequence."""
    used = [int(row["seed"]) for row in log if row.get("seed", "").isdigit()]
    return max(used) + 1 if used else seed_start


def write_backbone(out_dir: Path, name: str, pdb_text: str) -> Path:
    """Relabel chains, then write.

    The relabelling is the one piece of real work on this side: RFdiffusion
    returns a single renumbered chain, and every interface number computed
    later depends on chain A/C being the target and B the binder.
    """
    source = PDBFile.read(io.StringIO(pdb_text))
    array = source.get_structure(model=1)
    relabelled = nim.relabel_chains(array, CONTIGS)

    destination = out_dir / "backbones" / name
    out = PDBFile()
    out.set_structure(relabelled)
    out.write(str(destination))
    return destination


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default=DEFAULT_TARGET)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--key-file", type=Path, default=DEFAULT_KEY_FILE)
    parser.add_argument("--n", type=int, default=8, help="designs to generate")
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--diffusion-steps", type=int, default=DIFFUSION_STEPS)
    parser.add_argument(
        "--until",
        help="ISO UTC cutoff, e.g. 2026-11-09T22:00:00; overrides --n",
    )
    parser.add_argument(
        "--set-key",
        action="store_true",
        help="prompt for the API key, save it, and exit without generating",
    )
    args = parser.parse_args()

    if args.set_key:
        nim.prompt_for_key(args.key_file.expanduser())
        return 0

    target_path = RELAXED / args.target
    if not target_path.is_file():
        raise SystemExit(f"target not found: {target_path}")

    key = nim.read_api_key(args.key_file)
    out_dir = args.out_dir.expanduser()
    (out_dir / "backbones").mkdir(parents=True, exist_ok=True)

    cutoff = None
    if args.until:
        cutoff = datetime.fromisoformat(args.until).replace(tzinfo=UTC)

    log = read_log(out_dir)
    seed = next_seed(log, args.seed_start)
    batch = (max(int(r["batch"]) for r in log) + 1) if log else 0

    print(f"target  : {target_path.name}")
    print(f"contigs : {CONTIGS}")
    print(f"hotspots: {HOTSPOTS}   steps: {args.diffusion_steps}")
    print(f"noise   : {NOISE_SCALE}  <-- not settable via the API")
    print(f"out     : {out_dir}")
    print(f"starting at batch {batch}, seed {seed}")

    pdb_text = target_path.read_text()
    written = 0
    failures = 0

    while True:
        if cutoff is not None:
            remaining = (cutoff - datetime.now(UTC)).total_seconds()
            if remaining <= 0:
                print("cutoff reached")
                break
            print(f"\n[{remaining / 60:.0f} min left] design {written}, seed {seed}")
        else:
            if written >= args.n:
                break
            print(f"\ndesign {written + 1}/{args.n}, seed {seed}")

        payload = nim.build_payload(
            pdb_text, CONTIGS, HOTSPOTS, args.diffusion_steps, seed
        )
        try:
            design = nim.generate_one(key, payload)
        except nim.NimError as error:
            # One bad design should not end a production run, but a wrong
            # request would fail identically every time, so give up after a few.
            failures += 1
            print(f"  FAILED: {error}")
            if failures >= 3 and written == 0:
                raise SystemExit(
                    "three failures and nothing generated - the request itself "
                    "is being rejected, not the service being flaky"
                ) from error
            seed += 1
            continue

        name = f"nim{batch}_{written}.pdb"
        path = write_backbone(out_dir, name, design.pdb_text)
        append_log(out_dir, {
            "batch": batch,
            "seed": seed,
            "design": written,
            "filename": name,
            "contigs": CONTIGS,
            "hotspots": ",".join(HOTSPOTS),
            "noise_scale": NOISE_SCALE,
            "diffuser_T": args.diffusion_steps,
            "checkpoint": CHECKPOINT,
            "elapsed_s": round(design.elapsed_s, 1),
            "written_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        })
        server = f", server {design.server_elapsed_ms} ms" if design.server_elapsed_ms else ""
        print(f"  wrote {path.name} in {design.elapsed_s:.0f}s{server}")

        written += 1
        seed += 1

    print(f"\n{written} backbone(s) in {out_dir / 'backbones'}  ({failures} failure(s))")
    if written:
        print("Score them with:")
        print(f'  python scripts/05_ingest_backbones.py --drive-dir "{out_dir}"')
    return 0 if written else 1


if __name__ == "__main__":
    raise SystemExit(main())
