from collections import defaultdict

import structlog

from ref_builder.models import Molecule
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.schema import Segment, OTUSchema
from ref_builder.utils import IsolateName, Accession, IsolateNameType

logger = structlog.get_logger("otu.utils")


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
            versioned_accession = Accession.create_from_string(record.accession_version)
            isolates[isolate_name][versioned_accession] = record

    return isolates


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
            segment_name = record.source.organism

        return [Segment(name=segment_name, required=True, length=len(record.sequence))]

    segments = []
    for record in sorted(records, key=lambda record: record.accession):
        if record.source.segment:
            segments.append(
                Segment(
                    name=record.source.segment,
                    required=True,
                    length=len(record.sequence),
                ),
            )
        else:
            raise ValueError("No segment name found for multipartite OTU segment.")

    return segments


