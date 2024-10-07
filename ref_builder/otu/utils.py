import re
from collections import defaultdict
from enum import StrEnum
from uuid import UUID, uuid4

import structlog

from ref_builder.models import Molecule
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.plan import (
    MonopartitePlan,
    MultipartitePlan,
    SegmentPlan,
    SegmentRule,
    get_multipartite_segment_name,
)
from ref_builder.utils import Accession, IsolateName, IsolateNameType

logger = structlog.get_logger("otu.utils")


class DeleteRationale(StrEnum):
    """Default strings delineating reasons for resource deletion."""

    USER = "Requested by user"
    REFSEQ = "Superceded by RefSeq"


class RefSeqConflictError(ValueError):
    """Raised when a potential RefSeq replacement is found."""

    def __init__(
        self,
        message: str,
        isolate_id: UUID,
        isolate_name: IsolateName,
        accessions: list[str],
    ):
        super().__init__(message)

        self.isolate_id = isolate_id

        self.isolate_name = isolate_name

        self.accessions = accessions


def create_isolate_plan_from_records(
    records: list[NCBIGenbank],
    segments: list[SegmentPlan] | None = None,
) -> MonopartitePlan | MultipartitePlan | None:
    """Return a plan from a list of records representing an isolate."""
    if len(records) == 1:
        return MonopartitePlan(id=uuid4(), length=len(records[0].sequence))

    binned_records = group_genbank_records_by_isolate(records)
    if len(binned_records) > 1:
        logger.fatal(
            "More than one isolate found. Cannot create schema automatically.",
            bins=binned_records,
        )
        return None

    if segments is not None:
        return MultipartitePlan(id=uuid4(), segments=segments)

    segments = []
    for record in sorted(records, key=lambda record: record.accession):
        if not record.source.segment:
            raise ValueError("No segment name found for multipartite OTU segment.")

        segments.append(
            SegmentPlan(
                id=uuid4(),
                name=get_multipartite_segment_name(record),
                required=SegmentRule.REQUIRED,
                length=len(record.sequence),
            ),
        )

    if segments:
        return MultipartitePlan(id=uuid4(), segments=segments)

    return None


def group_genbank_records_by_isolate(
    records: list[NCBIGenbank],
) -> dict[IsolateName, dict[Accession, NCBIGenbank]]:
    """Indexes Genbank records by isolate name."""
    isolates = defaultdict(dict)

    for record in records:
        try:
            isolate_name = _get_isolate_name(record)

            if isolate_name is None:
                logger.debug(
                    "RefSeq record does not contain sufficient source data "
                    " for automatic inclusion. Add this record manually.",
                    accession=record.accession,
                    definition=record.definition,
                    source_data=record.source,
                )
                continue

            versioned_accession = Accession.from_string(record.accession_version)
            isolates[isolate_name][versioned_accession] = record

        except ValueError:
            continue

    return isolates


def parse_refseq_comment(comment: str) -> tuple[str, str]:
    """Parse a standard RefSeq comment."""
    if not comment:
        raise ValueError("Empty comment")

    refseq_pattern = re.compile(r"^(\w+ REFSEQ): [\w ]+. [\w ]+ ([\w]+).")

    match = refseq_pattern.search(comment)

    return match.group(1), match.group(2)


def get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    """Return relevant molecule metadata from one or more records."""
    if not records:
        raise IndexError("No records given")

    rep_record = None
    for record in records:
        if record.refseq:
            rep_record = record
            break
    if rep_record is None:
        rep_record = records[0]

    return Molecule.model_validate(
        {
            "strandedness": rep_record.strandedness.value,
            "type": rep_record.moltype.value,
            "topology": rep_record.topology.value,
        },
    )


def _get_isolate_name(record: NCBIGenbank) -> IsolateName | None:
    """Get the isolate name from a Genbank record."""
    if record.source.model_fields_set.intersection(
        {IsolateNameType.ISOLATE, IsolateNameType.STRAIN, IsolateNameType.CLONE},
    ):
        for source_type in IsolateNameType:
            if source_type in record.source.model_fields_set:
                return IsolateName(
                    type=IsolateNameType(source_type),
                    value=record.source.model_dump()[source_type],
                )

    if record.refseq:
        return None

    raise ValueError("Record does not contain sufficient source data for inclusion.")
