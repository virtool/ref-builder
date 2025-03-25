from uuid import UUID

from pydantic import ValidationError
from structlog import get_logger

from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU
from ref_builder.otu.models import OTU


logger = get_logger("otu.validate")


def get_validated_otu(unvalidated_otu: RepoOTU) -> OTU:
    """Get a validated OTU from a RepoOTU instance."""
    return OTU(
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


def check_otu_for_validity(unvalidated_otu: RepoOTU) -> bool:
    """Assure that an OTU can pass the validation standard."""
    try:
        get_validated_otu(unvalidated_otu)
    except ValidationError as exc:
        logger.error("OTU does not pass validation.", errors=exc.errors())
        return False

    return True


def check_repo_for_invalid_otus(repo: Repo) -> set[UUID]:
    """Run validation on all OTUs and return the Ids of bad OTUs."""
    invalid_otu_ids = set()

    for unvalidated_otu in repo.iter_otus():
        otu_logger = logger.bind(
            otu_id=str(unvalidated_otu.id),
            name=unvalidated_otu.name,
            taxid=unvalidated_otu.taxid,
        )
        if check_otu_for_validity(unvalidated_otu):
            otu_logger.debug(
                "OTU passed validation",
            )
        else:
            otu_logger.debug(
                "OTU did not pass validation",
            )
            invalid_otu_ids.add(unvalidated_otu.id)

    return invalid_otu_ids
