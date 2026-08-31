"""Immutable identity for authenticated Figure-6 claim artifacts."""

from typing import Final

from tools.measurement_claim_identity import EXTRACTION_COMMIT as EXTRACTION_COMMIT

ARTIFACT_KIND: Final = "benchmark-measurement-claims/v1"
INPUT_VERSION: Final = "benchmark-measurement-claim-input/v1"
RESULT_VERSION: Final = "benchmark-measurement-claim-result/v1"
MARKDOWN_NAME: Final = "figure6.md"
RAW_JSON_NAME: Final = "figure6.json"
PAYLOAD_NAME: Final = "benchmark-claims.json"
HISTORY_COMMIT: Final = "70049dd991e7d5c7b93d512785393a49a7f03564"
MC013_HISTORY_COMMIT: Final = "1e7d48e3d6b115d1ab4cb61c43ca21d2bc9bb6fd"
