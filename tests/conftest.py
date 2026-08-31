from __future__ import annotations

from pathlib import Path

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption("--held-publisher-pdf", type=Path)


@pytest.fixture
def held_publisher_pdf(request: pytest.FixtureRequest) -> Path:
    supplied = request.config.getoption("--held-publisher-pdf")
    if supplied is None:
        raise pytest.UsageError("--held-publisher-pdf is required by the live Held reference oracle")
    if not isinstance(supplied, Path):
        raise pytest.UsageError("--held-publisher-pdf must resolve to a filesystem path")
    return supplied.resolve(strict=True)
