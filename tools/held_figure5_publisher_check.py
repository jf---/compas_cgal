"""Check tracked Figure 5 publisher evidence against an official SVG."""

from __future__ import annotations

import argparse
from pathlib import Path

from benchmarks.held_figure5_publisher import InvalidFigure5PublisherEvidenceError
from benchmarks.held_figure5_publisher import extract_figure5_publisher_path
from benchmarks.held_figure5_publisher import load_figure5_publisher_path
from benchmarks.held_reference_cases import load_held_reference_case


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("svg", type=Path)
    args = parser.parse_args()
    if extract_figure5_publisher_path(args.svg, load_held_reference_case("figure5")) != load_figure5_publisher_path():
        raise InvalidFigure5PublisherEvidenceError("Tracked Figure 5 publisher evidence differs from the official SVG.")


if __name__ == "__main__":
    main()
