import warnings
from typing import Any

from pydantic import ValidationError
from structlog import get_logger

from ref_builder.otu.models import OTU
from ref_builder.resources import RepoOTU

logger = get_logger("otu.validate")


def check_otu_is_valid(unvalidated_otu: RepoOTU) -> bool:
    """Assure that an OTU can pass the validation standard."""
    try:
        with warnings.catch_warnings(record=True) as warning_list:
            validated_otu = OTU.model_validate(
                {
                    "id": unvalidated_otu.id,
                    "legacy_id": unvalidated_otu.legacy_id,
                    "name": unvalidated_otu.name,
                    "taxid": unvalidated_otu.taxid,
                    "acronym": unvalidated_otu.acronym,
                    "molecule": unvalidated_otu.molecule,
                    "plan": unvalidated_otu.plan,
                    "excluded_accessions": unvalidated_otu.excluded_accessions,
                    "isolates": unvalidated_otu.isolates,
                    "representative_isolate": unvalidated_otu.representative_isolate,
                }
            )

    except ValidationError as exc:
        logger.warning(
            "OTU does not pass validation.",
            id=str(unvalidated_otu.id),
            name=unvalidated_otu.name,
            taxid=unvalidated_otu.taxid,
        )

        otu_errors = otu_error_handler(exc)

        for error in otu_errors:
            logger.warning(
                "ValidationError",
                msg=error["msg"],
                type=error["type"],
                loc=error["loc"],
                input=error["input"] if not isinstance(error["input"], dict) else {},
            )
        return False

    if warning_list:
        logger.warning("Outstanding warnings found.", warning_count=len(warning_list))

        for warning_msg in warning_list:
            logger.warning(
                warning_msg.message,
                otu_id=str(validated_otu.id),
                warning_category=warning_msg.category.__name__,
            )

    return True


def otu_error_handler(e: ValidationError) -> list[dict[str, Any]]:
    """Add contextual metadata to OTU validation errors."""
    new_errors: list[dict[str, Any]] = e.errors()

    for error in new_errors:
        if not error["loc"]:
            if error["type"] == "segment_not_found":
                error["loc"] = (f"isolates[{error['ctx'].get('isolate_id')}]",)

            if error["type"] in ("sequence_too_short", "sequence_too_long"):
                error["loc"] = (
                    (
                        f"isolates[{error['ctx'].get('isolate_id')}]",
                        f"sequences[{error['ctx'].get('sequence_id', '')}]",
                    ),
                )

    return new_errors
