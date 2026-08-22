"""Run a named corpus, or the Figure 6 reproduction, and write its artifacts.

    python -m benchmarks.cli corpus --name all
    python -m benchmarks.cli corpus --name external --external-dir /path/to/profiles
    python -m benchmarks.cli figure6
    python -m benchmarks.cli figures

The corpus and Figure 6 commands produce different artifacts and are deliberately not one flag:
a corpus run emits `MeasurementRecord` rows over many instances, while the
Figure 6 reproduction emits a two-curve comparison over one pocket. Folding them
into a single `--corpus figure6` would put two unrelated shapes behind one name.
`figures` is a third shape again: it regenerates the published drawings of the
paths themselves, and exists so that no figure in the docs comes from a script
that is not in the repository.

THEY ALSO DEFAULT TO DIFFERENT PLACES, on purpose. A corpus report is dominated
by wall times, so it is a measurement OF A MACHINE and belongs in `build/`, never
committed. The Figure 6 comparison is lengths and angles, which are properties of
the geometry and reproduce anywhere, so it is a publishable artifact and defaults
into `docs/`.

Every sweep coordinate this module chooses is a named constant. The corpora are
meant to be re-run and compared across branches, so a number changed in passing
silently invalidates every earlier artifact it is compared against.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List
from typing import Optional
from typing import Sequence
from typing import Tuple

from benchmarks.degeneracy import degeneracy_corpus
from benchmarks.errors import BenchmarkError
from benchmarks.errors import MissingExternalDirectoryError
from benchmarks.errors import UnknownCorpusError
from benchmarks.external import RejectedProfile
from benchmarks.external import load_profiles
from benchmarks.families.analytic import arc_channel
from benchmarks.families.analytic import disk
from benchmarks.families.analytic import rectangle
from benchmarks.families.analytic import stadium
from benchmarks.families.complexity import arc_fraction_sweep
from benchmarks.families.complexity import ngon_sweep
from benchmarks.families.necks import PINCH_SWEEP_DEFAULT
from benchmarks.families.necks import pinch_sweep
from benchmarks.families.precision import DIGIT_SWEEP_DEFAULT
from benchmarks.families.precision import SCALE_SWEEP_DEFAULT
from benchmarks.families.precision import digit_sweep
from benchmarks.families.precision import scale_sweep
from benchmarks.families.topology import island_sweep
from benchmarks.figure6 import FIGURE6_CAPS
from benchmarks.figure6 import reference_pocket
from benchmarks.figure6 import run_figure6
from benchmarks.figure6 import write_figure6
from benchmarks.figures import DEFAULT_FIGURE_FORMATS
from benchmarks.figures import DEFAULT_FIGURES_OUT
from benchmarks.figures import regenerate_every_figure
from benchmarks.mathsm import SPACING_SWEEP_TOOL_DIAMETERS
from benchmarks.palette import Theme
from benchmarks.report import write_report
from benchmarks.runner import run_corpus
from benchmarks.spec import PocketSpec

CORPUS_NAMES: Tuple[str, ...] = ("smoke", "analytic", "complexity", "necks", "precision", "topology", "degeneracy", "external", "all")

# What `all` expands to. `smoke` is a subset of `analytic` and `external` needs an
# operator-supplied directory, so neither belongs in the aggregate.
AGGREGATE_CORPUS_NAMES: Tuple[str, ...] = ("analytic", "complexity", "necks", "precision", "topology", "degeneracy")

# Corpus reports carry wall times and are therefore machine-specific; they stay
# out of the published site. The Figure 6 comparison is geometry and reproduces
# anywhere, so it is written where the docs can link to it.
DEFAULT_CORPUS_OUT = Path("build/benchmarks")
DEFAULT_FIGURE6_OUT = Path("docs/benchmarks")

DEFAULT_TOOL_DIAMETER = 1.0
DEFAULT_CAP_DEG = 120.0

# --- smoke: one instance, sized so a CI run is a second and not a minute -------
SMOKE_WIDTH = 8.0
SMOKE_HEIGHT = 6.0
SMOKE_TOOL_DIAMETER = 2.0

# --- analytic: four shapes whose clearance is derivable in closed form ---------
ANALYTIC_DISK_RADIUS = 8.0
ANALYTIC_RECT_WIDTH = 20.0
ANALYTIC_RECT_HEIGHT = 10.0
ANALYTIC_STADIUM_LENGTH = 20.0
ANALYTIC_STADIUM_HALF_WIDTH = 3.0
ANALYTIC_CHANNEL_GUIDE_RADIUS = 10.0
ANALYTIC_CHANNEL_HALF_WIDTH = 2.0
ANALYTIC_CHANNEL_SWEEP_DEG = 90.0

# --- complexity: boundary element count, then boundary element KIND ------------
# The k sweep holds area fixed so only the element count moves; the arc sweep
# holds the count fixed so only the fraction of circular sides moves.
NGON_SIDE_COUNTS: Tuple[int, ...] = (3, 6, 12, 32, 128)
NGON_AREA = 100.0
ARC_FRACTION_SIDES = 16
ARC_FRACTION_RATIOS: Tuple[float, ...] = (0.0, 0.25, 0.5, 1.0)
ARC_FRACTION_RADIUS = 10.0

# --- precision: coordinate magnitude, then significant decimals ----------------
PRECISION_SIDE_COUNT = 12
PRECISION_SEED = 1
PRECISION_SCALE_DECIMALS = 6
PRECISION_DIGIT_RADIUS = 10.0

# --- topology: island count at a fixed outer boundary -------------------------
ISLAND_GRID_COUNTS: Tuple[Tuple[int, int], ...] = ((1, 1), (2, 2), (3, 3), (4, 4))

EXIT_OK = 0
EXIT_CORPUS_ERROR = 2


@dataclass(frozen=True)
class CorpusBuild:
    """The instances a corpus name resolved to, and what it could not use.

    Attributes:
        name: The corpus name that was requested.
        specs: Instances to measure, in corpus order.
        rejected: Profiles an external corpus could not turn into instances.
            Always empty for an authored corpus, and never dropped for a loaded
            one -- a thinner corpus that reports as a full one is the failure
            mode this field exists to prevent.
    """

    name: str
    specs: Tuple[PocketSpec, ...]
    rejected: Tuple[RejectedProfile, ...]


def build_corpus(name: str, tool_diameter: float, tea_cap_deg: float, external_dir: Optional[Path] = None) -> CorpusBuild:
    """Assemble the named corpus.

    Args:
        name: One of `CORPUS_NAMES`.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.
        external_dir: Directory of third-party profiles, for the external corpus.

    Returns:
        The build.

    Raises:
        UnknownCorpusError: The corpus name is not in `CORPUS_NAMES`.
        MissingExternalDirectoryError: `external` was requested without a
            directory.
        ExternalCorpusError: An external directory is missing or malformed.
    """
    if name == "external":
        if external_dir is None:
            raise MissingExternalDirectoryError("The external corpus needs --external-dir; nothing in this repository supplies third-party geometry.")
        corpus = load_profiles(external_dir, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg)
        return CorpusBuild(name=name, specs=corpus.specs, rejected=corpus.rejected)
    if name == "all":
        specs: List[PocketSpec] = []
        for sub in AGGREGATE_CORPUS_NAMES:
            specs.extend(build_corpus(sub, tool_diameter, tea_cap_deg).specs)
        return CorpusBuild(name=name, specs=tuple(specs), rejected=())
    return CorpusBuild(name=name, specs=tuple(_authored_specs(name, tool_diameter, tea_cap_deg)), rejected=())


def _authored_specs(name: str, tool_diameter: float, tea_cap_deg: float) -> List[PocketSpec]:
    """The instances of one authored corpus.

    Args:
        name: A corpus name other than `external` and `all`.
        tool_diameter: Cutter diameter.
        tea_cap_deg: Engagement cap in degrees.

    Returns:
        The instances, in corpus order.

    Raises:
        UnknownCorpusError: The corpus name is not in `CORPUS_NAMES`.
    """
    if name == "smoke":
        return [rectangle(width=SMOKE_WIDTH, height=SMOKE_HEIGHT, tool_diameter=SMOKE_TOOL_DIAMETER, tea_cap_deg=tea_cap_deg)]
    if name == "analytic":
        return [
            disk(radius=ANALYTIC_DISK_RADIUS, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
            rectangle(width=ANALYTIC_RECT_WIDTH, height=ANALYTIC_RECT_HEIGHT, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
            stadium(straight_length=ANALYTIC_STADIUM_LENGTH, half_width=ANALYTIC_STADIUM_HALF_WIDTH, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg),
            arc_channel(
                guide_radius=ANALYTIC_CHANNEL_GUIDE_RADIUS,
                half_width=ANALYTIC_CHANNEL_HALF_WIDTH,
                sweep_deg=ANALYTIC_CHANNEL_SWEEP_DEG,
                tool_diameter=tool_diameter,
                tea_cap_deg=tea_cap_deg,
            ),
        ]
    if name == "complexity":
        return ngon_sweep(ks=NGON_SIDE_COUNTS, area=NGON_AREA, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg) + arc_fraction_sweep(
            n=ARC_FRACTION_SIDES, ratios=ARC_FRACTION_RATIOS, radius=ARC_FRACTION_RADIUS, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg
        )
    if name == "necks":
        return pinch_sweep(pinches=PINCH_SWEEP_DEFAULT, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg)
    if name == "precision":
        return scale_sweep(k=PRECISION_SIDE_COUNT, scales=SCALE_SWEEP_DEFAULT, decimals=PRECISION_SCALE_DECIMALS, seed=PRECISION_SEED, tea_cap_deg=tea_cap_deg) + digit_sweep(
            k=PRECISION_SIDE_COUNT,
            radius=PRECISION_DIGIT_RADIUS,
            decimal_counts=DIGIT_SWEEP_DEFAULT,
            seed=PRECISION_SEED,
            tool_diameter=tool_diameter,
            tea_cap_deg=tea_cap_deg,
        )
    if name == "topology":
        return island_sweep(counts=ISLAND_GRID_COUNTS, tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg)
    if name == "degeneracy":
        return degeneracy_corpus(tool_diameter=tool_diameter, tea_cap_deg=tea_cap_deg)
    raise UnknownCorpusError(f"Unknown corpus {name!r}; expected one of {CORPUS_NAMES}.")


def _run_corpus_command(args: argparse.Namespace) -> int:
    """Measure a corpus and write its report.

    Args:
        args: Parsed `corpus` arguments.

    Returns:
        `EXIT_OK`.
    """
    build = build_corpus(args.name, args.tool_diameter, args.cap_deg, args.external_dir)
    for rejection in build.rejected:
        print(f"rejected {rejection.name}: {rejection.reason}", file=sys.stderr)
    records = run_corpus(build.specs, collect_digits=not args.no_digits)
    md_path, json_path = write_report(records, args.out)
    failed = sum(1 for r in records if r.error is not None)
    print(f"corpus {build.name}: {len(records)} measured ({failed} failed), {len(build.rejected)} rejected before measuring")
    print(f"wrote {md_path} and {json_path}")
    return EXIT_OK


def _run_figure6_command(args: argparse.Namespace) -> int:
    """Reproduce Figure 6 and write its artifacts.

    Args:
        args: Parsed `figure6` arguments.

    Returns:
        `EXIT_OK`.
    """
    run = run_figure6(reference_pocket(), caps=tuple(args.caps), spacings_tool_diameters=tuple(args.spacings))
    md_path, json_path = write_figure6(run, args.out)
    print(f"figure6 on {run.spec.name}: {len(run.points)} caps against {len(run.trials)} constant-spacing trials")
    print(f"wrote {md_path} and {json_path}")
    return EXIT_OK


def _run_figures_command(args: argparse.Namespace) -> int:
    """Redraw every published figure -- the tool-path comparison and the quality set.

    Args:
        args: Parsed `figures` arguments.

    Returns:
        `EXIT_OK`.
    """
    written = regenerate_every_figure(args.out, formats=tuple(args.formats), themes=(Theme(args.theme),))
    for path in written:
        print(f"wrote {path}")
    return EXIT_OK


def build_parser() -> argparse.ArgumentParser:
    """The command-line interface.

    Returns:
        The parser.
    """
    parser = argparse.ArgumentParser(prog="benchmarks.cli", description="Run a pocket-machining benchmark corpus, or reproduce Held & Pfeiffer's Figure 6.")
    commands = parser.add_subparsers(dest="command", required=True)

    corpus = commands.add_parser("corpus", help="Measure a named corpus and write its report.")
    corpus.add_argument("--name", choices=CORPUS_NAMES, default="smoke", help="Which corpus to measure.")
    corpus.add_argument("--tool-diameter", type=float, default=DEFAULT_TOOL_DIAMETER, help="Cutter diameter.")
    corpus.add_argument("--cap-deg", type=float, default=DEFAULT_CAP_DEG, help="Engagement cap in degrees.")
    corpus.add_argument("--out", type=Path, default=DEFAULT_CORPUS_OUT, help="Directory for the report artifacts.")
    corpus.add_argument("--external-dir", type=Path, default=None, help="Directory of operator-prepared third-party profiles.")
    corpus.add_argument("--no-digits", action="store_true", help="Skip the untimed exact-coordinate diagnostic pass.")
    corpus.set_defaults(handler=_run_corpus_command)

    figure6 = commands.add_parser("figure6", help="Reproduce Figure 6 on the reference pocket (minutes, not seconds).")
    figure6.add_argument("--out", type=Path, default=DEFAULT_FIGURE6_OUT, help="Directory for the figure artifacts.")
    figure6.add_argument("--caps", type=float, nargs="+", default=list(FIGURE6_CAPS), help="Engagement caps in degrees, in plotting order.")
    figure6.add_argument("--spacings", type=float, nargs="+", default=list(SPACING_SWEEP_TOOL_DIAMETERS), help="Baseline trial spacings, in tool diameters.")
    figure6.set_defaults(handler=_run_figure6_command)

    figures = commands.add_parser("figures", help="Redraw every published tool-path figure (seconds: it generates the paths, it does not audit them).")
    figures.add_argument("--out", type=Path, default=DEFAULT_FIGURES_OUT, help="Directory for the figure files.")
    figures.add_argument("--formats", nargs="+", default=list(DEFAULT_FIGURE_FORMATS), help="File formats to write, one file each.")
    figures.add_argument("--theme", choices=[theme.value for theme in Theme], default=Theme.LIGHT.value, help="Which surface the figure is drawn for.")
    figures.set_defaults(handler=_run_figures_command)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run one command.

    A corpus failure is reported and exits non-zero rather than raising: the
    caller asked for something this repository cannot supply (an absent external
    directory, an unknown name), which is a usage error, not a defect. Anything
    else propagates, because a defect must not be turned into an exit code.

    Args:
        argv: Command-line arguments; defaults to `sys.argv[1:]`.

    Returns:
        `EXIT_OK`, or `EXIT_CORPUS_ERROR` when the corpus could not be assembled.
    """
    args = build_parser().parse_args(argv)
    try:
        exit_code: int = args.handler(args)
    except BenchmarkError as exc:
        print(f"error: {type(exc).__name__}: {exc}", file=sys.stderr)
        return EXIT_CORPUS_ERROR
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
