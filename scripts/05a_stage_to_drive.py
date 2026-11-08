#!/usr/bin/env python
"""Stage 05a - copy design targets into a Google Drive folder for Colab.

RFdiffusion needs a GPU this project does not have locally, so stage 05 runs in
Colab. Google Drive for Desktop is mounted on this machine, which makes the
transfer a plain file copy rather than an API upload: no OAuth client, no
credentials, and nothing embedded in the notebook.

The notebook reads ``TARGETS.tsv`` written here and asserts the SHA-256 of the
file it loads against the value recorded at staging time. That is what ties the
generated backbones to a specific stage-04 output, so a silently wrong or
half-synced target fails loudly instead of producing plausible nonsense.

Usage:
    python scripts/05a_stage_to_drive.py --drive-dir "<...>/My Drive/kras_pmhc_stage05"
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RELAXED = ROOT / "data" / "targets_relaxed"

#: The single agreed Route A target. Copy L of 7OW3 has no backbone gap in the
#: cropped region, and chi1 -70 is both the most populated rotamer well and the
#: most exposed carboxylate. ``_crop`` is HLA alpha1/alpha2 (chain A, residues
#: 1-180) plus the peptide (chain C, 1-10) - the target RFdiffusion sees.
DEFAULT_TARGETS = ("7ow3_L_G6D_chi1m70_chi2p120_relaxed_crop.pdb",)

MANIFEST_NAME = "TARGETS.tsv"
MANIFEST_FIELDS = ("filename", "sha256", "n_bytes", "source_path", "staged_utc")

CHUNK = 1 << 20


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_targets(names: list[str]) -> list[Path]:
    """Accept either a bare filename in data/targets_relaxed/ or a full path."""
    paths = []
    for name in names:
        candidate = Path(name)
        path = candidate if candidate.is_absolute() or candidate.exists() else RELAXED / name
        if not path.is_file():
            raise SystemExit(f"target not found: {name} (looked at {path})")
        paths.append(path.resolve())
    return paths


def stage(source: Path, drive_dir: Path, force: bool) -> dict[str, str]:
    """Copy one target into the Drive folder and verify it arrived intact."""
    destination = drive_dir / source.name
    expected = sha256(source)

    if destination.exists() and not force:
        existing = sha256(destination)
        if existing == expected:
            print(f"  unchanged: {source.name}")
        else:
            raise SystemExit(
                f"{destination} already exists with different content "
                f"({existing[:16]} != {expected[:16]}); pass --force to overwrite"
            )
    else:
        shutil.copy2(source, destination)
        print(f"  copied:    {source.name}")

    # Verify after the copy, not before: Drive for Desktop writes through a
    # synced filesystem, and a truncated or failed write must not pass silently.
    if not destination.is_file():
        raise SystemExit(f"{destination} missing after copy")
    landed = sha256(destination)
    if landed != expected:
        raise SystemExit(
            f"{destination} checksum mismatch after copy: "
            f"{landed[:16]} != {expected[:16]}"
        )

    return {
        "filename": source.name,
        "sha256": expected,
        "n_bytes": str(source.stat().st_size),
        "source_path": str(source.relative_to(ROOT)),
        "staged_utc": datetime.now(UTC).isoformat(timespec="seconds"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--drive-dir",
        required=True,
        type=Path,
        help="destination folder inside the mounted Google Drive",
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=list(DEFAULT_TARGETS),
        help="target PDB filenames (in data/targets_relaxed/) or full paths",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite a destination file whose content differs",
    )
    args = parser.parse_args()

    sources = resolve_targets(args.targets)

    drive_dir: Path = args.drive_dir.expanduser()
    drive_dir.mkdir(parents=True, exist_ok=True)
    if not drive_dir.is_dir():
        raise SystemExit(f"not a directory: {drive_dir}")

    print(f"staging {len(sources)} target(s) -> {drive_dir}")
    rows = [stage(source, drive_dir, args.force) for source in sources]

    manifest = drive_dir / MANIFEST_NAME
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  wrote:     {MANIFEST_NAME}")

    print()
    for row in rows:
        print(f"{row['filename']}  sha256 {row['sha256'][:16]}...  {row['n_bytes']} bytes")

    print(
        "\nWait for Google Drive to finish syncing before running the notebook.\n"
        "The notebook re-checks these checksums and will fail if a file is\n"
        "missing or partially uploaded."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
