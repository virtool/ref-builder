from pydantic import ValidationError
from structlog import get_logger

from ref_builder.resources import RepoOTU
from ref_builder.otu.models import OTU


logger = get_logger("otu.validate")


def validate_otu(unvalidated_otu: RepoOTU) -> bool:
    """Assure that an OTU can pass the validation standard."""
    try:
        validated_otu = OTU(
            id=unvalidated_otu.id,
            legacy_id=unvalidated_otu.legacy_id,
            name=unvalidated_otu.name,
            taxid=unvalidated_otu.taxid,
            acronym=unvalidated_otu.acronym,
            molecule=unvalidated_otu.molecule,
            plan=unvalidated_otu.plan,
            isolates=unvalidated_otu.isolates,
            excluded_accessions=unvalidated_otu.excluded_accessions,
            representative_isolate=unvalidated_otu.representative_isolate,
        )
    except ValidationError as exc:
        logger.warning("OTU does not pass validation.", errors=exc.errors())
        return False

    return True
