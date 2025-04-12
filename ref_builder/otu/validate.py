from typing import Any

from pydantic import ValidationError
from structlog import get_logger

from ref_builder.otu.models import OTU
from ref_builder.resources import RepoOTU

logger = get_logger("otu.validate")


def check_otu_is_valid(unvalidated_otu: RepoOTU) -> bool:
    """Assure that an OTU can pass the validation standard."""
    try:
        OTU(
            id=unvalidated_otu.id,
            legacy_id=unvalidated_otu.legacy_id,
            name=unvalidated_otu.name,
            taxid=unvalidated_otu.taxid,
            acronym=unvalidated_otu.acronym,
            molecule=unvalidated_otu.molecule,
            plan=unvalidated_otu.plan,
            excluded_accessions=unvalidated_otu.excluded_accessions,
            representative_isolate=unvalidated_otu.representative_isolate,
            isolates=unvalidated_otu.isolates,
        )

    except ValidationError as exc:
        logger.warning(
            "OTU does not pass validation.",
            id=unvalidated_otu.id,
            name=unvalidated_otu.name,
            taxid=unvalidated_otu.taxid,
        )

        otu_errors = otu_error_handler(exc)

        for error in otu_errors:
            logger.warning(
                "ValidationError: " + error["msg"],
                type=error["type"],
                loc=error["loc"],
                input=error["input"] if not isinstance(error["input"], dict) else {},
            )
        return False

    return True


def otu_error_handler(e: ValidationError):
    new_errors: list[dict[str, Any]] = e.errors()

    for error in new_errors:
        if not error["loc"]:
            if error["type"] in ("sequence_too_short", "sequence_too_long"):
                error["loc"] = (
                    f"isolates[{error['ctx'].get('isolate_id')}]",
                    f"sequences[{error['ctx'].get('sequence_id', '')}]"),

    return new_errors
