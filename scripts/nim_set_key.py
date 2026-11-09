#!/usr/bin/env python
"""Store the NVIDIA NIM API key, asked for at a hidden prompt.

Deliberately takes no arguments: it is run by hand in a terminal, and an
argument parser is one more thing to get wrong at a keyboard. A key given as a
command-line argument would also land in the shell history and the process
table, which is the thing this script exists to avoid.

    python scripts/nim_set_key.py

The key is written to .nim_api_key with mode 0600 and is gitignored.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kras_pmhc import nim  # noqa: E402


def main() -> int:
    try:
        nim.prompt_for_key(ROOT / ".nim_api_key")
    except nim.NimError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
