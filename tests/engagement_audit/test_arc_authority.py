"""Single-path guard for exact audit arc depletion."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_authoritative_arc_sources_never_call_legacy_sweep() -> None:
    authoritative = (
        ROOT / "src" / "audit_arc_motion_2.cpp",
        ROOT / "src" / "exact_circle_chart_2.cpp",
        ROOT / "src" / "exact_depletion_2.cpp",
    )

    for source in authoritative:
        assert "subtract_arc_sweep" not in source.read_text()


def test_continuous_tea_uses_shared_frozen_circle_atlas() -> None:
    source = (ROOT / "src" / "continuous_tea_2" / "parameter_charts.cpp").read_text()

    assert '#include "../exact_circle_chart_atlas_2.h"' in source
    assert "exact_circle_chart_record" in source
    assert "quarter_circle_numerators" not in source


def test_audit_identity_reuses_authoritative_ccan_rational_encoder() -> None:
    source = (ROOT / "src" / "audit_digest_2.h").read_text()

    assert "canonical_encode_rational" in source
    assert "CCAN" not in source


def test_frozen_atlas_is_the_only_quarter_coefficient_table() -> None:
    atlas = ROOT / "src" / "exact_circle_chart_atlas_2.h"
    consumers = (
        ROOT / "src" / "exact_circle_chart_2.cpp",
        ROOT / "src" / "continuous_tea_2" / "parameter_charts.cpp",
    )

    assert "{1, 0, -1}" in atlas.read_text()
    for consumer in consumers:
        text = consumer.read_text()
        assert "{1, 0, -1}" not in text
        assert "EXACT_CIRCLE_CHART_ATLAS" in text or "exact_circle_chart_record" in text
