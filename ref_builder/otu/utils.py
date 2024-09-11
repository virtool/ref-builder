import re
from collections import defaultdict
from enum import StrEnum
from uuid import UUID, uuid4

import structlog

from ref_builder.models import Molecule
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.schema import OTUSchema, Segment
from ref_builder.utils import Accession, IsolateName, IsolateNameType

logger = structlog.get_logger("otu.utils")


class DeleteRationale(StrEnum):
    """Default strings delineating reasons for resource deletion."""
    USER = "Requested by user"
    REFSEQ = "Superceded by RefSeq"


class RefSeqConflictError(ValueError):
    """Raised when a potential RefSeq replacement is found."""
    def __init__(
        self, message, isolate_id: UUID, isolate_name: IsolateName, accessions: list[str]
    ):
        super().__init__(message)

        self.isolate_id = isolate_id

        self.isolate_name = isolate_name

        self.accessions = accessions


def create_schema_from_records(
    records: list[NCBIGenbank],
    segments: list[Segment] | None = None,
) -> OTUSchema | None:
    molecule = _get_molecule_from_records(records)

    binned_records = group_genbank_records_by_isolate(records)
    if len(binned_records) > 1:
        logger.fatal(
            "More than one isolate found. Cannot create schema automatically.",
            bins=binned_records,
        )
        return None

    if segments is None:
        segments = _get_segments_from_records(records)
        if segments:
            return OTUSchema(molecule=molecule, segments=segments)

    return None


def group_genbank_records_by_isolate(
    records: list[NCBIGenbank],
) -> dict[IsolateName, dict[Accession, NCBIGenbank]]:
    """Indexes Genbank records by isolate name"""
    isolates = defaultdict(dict)

    for record in records:
        if (isolate_name := _get_isolate_name(record)) is not None:
            versioned_accession = Accession.from_string(record.accession_version)
            isolates[isolate_name][versioned_accession] = record

    return isolates


def parse_refseq_comment(comment: str) -> tuple[str, str]:
    """Parse a standard RefSeq comment."""
    if not comment:
        raise ValueError("Empty comment")

    refseq_pattern = re.compile(r"^(\w+ REFSEQ): [\w ]+. [\w ]+ ([\w]+).")

    match = refseq_pattern.search(comment)

    return match.group(1), match.group(2)


def _get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    """Return relevant molecule metadata from one or more records"""
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
    """Get the isolate name from a Genbank record"""
    record_logger = logger.bind(
        accession=record.accession,
        definition=record.definition,
        source_data=record.source,
    )

    if record.source.model_fields_set.intersection(
        {IsolateNameType.ISOLATE, IsolateNameType.STRAIN, IsolateNameType.CLONE},
    ):
        for source_type in IsolateNameType:
            if source_type in record.source.model_fields_set:
                return IsolateName(
                    type=IsolateNameType(source_type),
                    value=record.source.model_dump()[source_type],
                )

    elif record.refseq:
        record_logger.debug(
            "RefSeq record does not contain sufficient source data. "
            + "Edit before inclusion.",
        )

        return IsolateName(
            type=IsolateNameType(IsolateNameType.REFSEQ),
            value=record.accession,
        )

    record_logger.debug("Record does not contain sufficient source data for inclusion.")

    return None


def _get_segments_from_records(records: list[NCBIGenbank]) -> list[Segment]:
    if len(records) == 1:
        record = records[0]

        if record.source.segment != "":
            segment_name = record.source.segment
        else:
            segment_name = record.source.mol_type

        segment_id = uuid4()
        return [
            Segment(
                id=segment_id, name=segment_name, required=True, length=len(record.sequence)
            )
        ]

    segments = []
    for record in sorted(records, key=lambda record: record.accession):
        if record.source.segment:
            segment_id = uuid4()
            segments.append(
                Segment(
                    id=segment_id,
                    name=record.source.segment,
                    required=True,
                    length=len(record.sequence),
                ),
            )
        else:
            raise ValueError("No segment name found for multipartite OTU segment.")

    return segments
