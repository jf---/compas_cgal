"""Machine and build identity shared by the performance records, and the cross-machine time policy."""

from __future__ import annotations

import pathlib
import platform
import subprocess
from dataclasses import dataclass
from datetime import datetime
from datetime import timezone

from benchmarks.errors import BenchmarkError
from compas_cgal import __version__ as PACKAGE_VERSION

REPOSITORY = pathlib.Path(__file__).resolve().parents[1]

# Run-to-run spread of one workload on a quiet machine was measured under
# 1.5x on 2026-09-08; 3x leaves room for a loaded box while catching the
# order-of-magnitude regressions this project has produced. Across machines
# only order of magnitude is comparable, so the factor is 10x and each record
# carries its own absolute ceilings.
SAME_MACHINE_FACTOR = 3.0
CROSS_MACHINE_FACTOR = 10.0


class BuildIdentityError(BenchmarkError):
    """The repository identity could not be read from git."""


@dataclass(frozen=True)
class MachineIdentity:
    machine: str
    system: str
    cpu_brand: str
    python: str

    def matches(self, other: MachineIdentity) -> bool:
        """Same hardware class: the same-machine factors apply."""
        return (self.machine, self.system, self.cpu_brand) == (other.machine, other.system, other.cpu_brand)


@dataclass(frozen=True)
class BuildIdentity:
    commit: str
    dirty: bool
    package_version: str
    recorded_at: str


def cpu_brand() -> str:
    """Marketing name of the CPU, or 'unknown' where the platform does not publish one."""
    system = platform.system()
    if system == "Darwin":
        return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, check=True).stdout.strip()
    if system == "Linux":
        for line in pathlib.Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
        return "unknown"
    return platform.processor() or "unknown"


def machine_identity() -> MachineIdentity:
    return MachineIdentity(machine=platform.machine(), system=platform.system(), cpu_brand=cpu_brand(), python=platform.python_version())


def build_identity(repository: pathlib.Path = REPOSITORY) -> BuildIdentity:
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repository, capture_output=True, text=True, check=True).stdout.strip()
        status = subprocess.run(["git", "status", "--porcelain"], cwd=repository, capture_output=True, text=True, check=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as error:
        raise BuildIdentityError(f"git identity unavailable in {repository}: {error}") from error
    return BuildIdentity(commit=commit, dirty=bool(status.strip()), package_version=PACKAGE_VERSION, recorded_at=datetime.now(timezone.utc).isoformat(timespec="seconds"))


def time_factor(current: MachineIdentity, baseline: MachineIdentity) -> tuple[float, str]:
    """The factor a current time may exceed the baseline by, and the regime it came from."""
    if current.matches(baseline):
        return SAME_MACHINE_FACTOR, "same machine"
    return CROSS_MACHINE_FACTOR, f"cross-machine ({current.cpu_brand} vs {baseline.cpu_brand})"
