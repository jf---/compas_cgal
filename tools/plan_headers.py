"""Invariant I4: every plan carries a status header.

Status lives in the header and is derived from artifacts; bodies are immutable
planning history. This lint only asserts the header exists. Truthfulness is
the audit's job, not grep's.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import cast

HEADER_PREFIX = "> **status:"
SEARCH_LINES = 12
DEFAULT_ALLOW_PREFIXES: Tuple[str, ...] = ("2026-08-23-auditor-convergence",)


def missing_headers(
    plans_dir: Path,
    allow_prefixes: Sequence[str] = (),
) -> List[Path]:
    """Return plans lacking a status header in their first 12 lines.

    Args:
        plans_dir: Directory containing plan Markdown files.
        allow_prefixes: Plan-name prefixes owned by another governance wave.

    Returns:
        Sorted paths of non-exempt plans without an early status header.
    """
    exempt_prefixes = tuple(allow_prefixes)
    missing: List[Path] = []
    for plan in sorted(plans_dir.glob("*.md")):
        if plan.name.startswith(exempt_prefixes):
            continue
        head = plan.read_text().splitlines()[:SEARCH_LINES]
        if not any(line.startswith(HEADER_PREFIX) for line in head):
            missing.append(plan)
    return missing


def main(argv: Optional[List[str]] = None) -> int:
    """Run the plan-header check."""
    parser = argparse.ArgumentParser(description="Every plan must carry a '> **status:' header.")
    parser.add_argument(
        "--plans-dir",
        type=Path,
        default=Path("docs/superpowers/plans"),
    )
    parser.add_argument(
        "--allow-prefix",
        action="append",
        default=None,
        help=("Additional plan-name prefix exempt until its owner stamps it; repeat for multiple prefixes."),
    )
    args = parser.parse_args(argv)
    cli_allow_prefixes = cast(Optional[List[str]], args.allow_prefix)
    explicit_prefixes: Tuple[str, ...] = () if cli_allow_prefixes is None else tuple(cli_allow_prefixes)
    allow_prefixes = DEFAULT_ALLOW_PREFIXES + explicit_prefixes
    missing = missing_headers(args.plans_dir, allow_prefixes=allow_prefixes)
    for plan in missing:
        print(f"missing status header: {plan}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
