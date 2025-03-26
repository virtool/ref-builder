from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.create import create_otu_with_taxid
from ref_builder.otu.promote import replace_accessions_from_records, promote_otu_accessions_from_records
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate
from tests.fixtures.factories import IsolateFactory


def test_multi_link(empty_repo):
    with empty_repo.lock():
        otu_init = create_otu_with_taxid(
            empty_repo, 2164102, ["MF062136", "MF062137", "MF062138"], acronym=""
        )

    assert otu_init.accessions == {"MF062136", "MF062137", "MF062138"}

    segment_l = otu_init.plan.get_segment_by_name_key("L")

    segment_l_rep_sequence = otu_init.get_sequence_by_accession("MF062136")

    otu_id = otu_init.id

    mock_isolate = RepoIsolate.model_validate(
        IsolateFactory.build_on_plan(otu_init.plan).model_dump()
    )

    for iterator in range(len(mock_isolate.sequences)):
        mock_sequence = mock_isolate.sequences[iterator]

        if mock_sequence.segment == segment_l.id:
            mock_isolate.sequences[iterator] = segment_l_rep_sequence.copy()

    mock_isolate_accessions = set(
        sequence.accession.key for sequence in mock_isolate.sequences
    )

    with empty_repo.lock(), empty_repo.use_transaction():
        isolate_init = empty_repo.create_isolate(
            otu_id, legacy_id=None, name=mock_isolate.name
        )

        empty_repo.link_sequence(
            otu_id,
            isolate_init.id,
            segment_l_rep_sequence.id,
        )

        for mock_sequence in mock_isolate.sequences:
            # skip segment L
            if mock_sequence.segment == segment_l.id:
                continue

            sequence_init = empty_repo.create_sequence(
                otu_id,
                accession=str(mock_sequence.accession),
                definition=mock_sequence.definition,
                legacy_id=None,
                segment=mock_sequence.segment,
                sequence=mock_sequence.sequence,
            )

            empty_repo.link_sequence(
                otu_id,
                isolate_init.id,
                sequence_init.id,
            )

    otu_after_mock_isolate = empty_repo.get_otu(otu_id)

    refseq_records = NCBIClient(ignore_cache=False).fetch_genbank_records(
        ["NC_055390"]
    )

    with empty_repo.lock():
        replace_accessions_from_records(
            empty_repo,
            otu_after_mock_isolate,
            record_by_replaceable_sequence_id={segment_l_rep_sequence.id: refseq_records[0]},
        )

def test_multi_linked_promotion(empty_repo: Repo):
    with empty_repo.lock():
        otu_init = create_otu_with_taxid(
            empty_repo, 2164102, ["MF062136", "MF062137", "MF062138"], acronym=""
        )

    assert otu_init.accessions == {"MF062136", "MF062137", "MF062138"}

    segment_l = otu_init.plan.get_segment_by_name_key("L")

    segment_l_rep_sequence = otu_init.get_sequence_by_accession("MF062136")

    otu_id = otu_init.id

    mock_isolate = RepoIsolate.model_validate(
        IsolateFactory.build_on_plan(otu_init.plan).model_dump()
    )

    for iterator in range(len(mock_isolate.sequences)):
        mock_sequence = mock_isolate.sequences[iterator]

        if mock_sequence.segment == segment_l.id:
            mock_isolate.sequences[iterator] = segment_l_rep_sequence.copy()

    mock_isolate_accessions = set(
        sequence.accession.key for sequence in mock_isolate.sequences
    )

    with empty_repo.lock(), empty_repo.use_transaction():
        isolate_init = empty_repo.create_isolate(
            otu_id, legacy_id=None, name=mock_isolate.name
        )

        empty_repo.link_sequence(
            otu_id,
            isolate_init.id,
            segment_l_rep_sequence.id,
        )

        for mock_sequence in mock_isolate.sequences:
            # skip segment L
            if mock_sequence.segment == segment_l.id:
                continue

            sequence_init = empty_repo.create_sequence(
                otu_id,
                accession=str(mock_sequence.accession),
                definition=mock_sequence.definition,
                legacy_id=None,
                segment=mock_sequence.segment,
                sequence=mock_sequence.sequence,
            )

            empty_repo.link_sequence(
                otu_id,
                isolate_init.id,
                sequence_init.id,
            )

    otu_after_mock_isolate = empty_repo.get_otu(otu_id)

    assert isolate_init.id in otu_after_mock_isolate.isolate_ids

    assert (
        otu_after_mock_isolate.accessions
        == {"MF062136", "MF062137", "MF062138"} | mock_isolate_accessions
    )

    refseq_records = NCBIClient(ignore_cache=False).fetch_genbank_records(
        ["NC_055390", "NC_055391", "NC_055392"],
    )

    new_refseq_isolate = update_isolate_from_records(
        repo=empty_repo,
        otu=otu_after_mock_isolate,
        isolate_id=otu_init.representative_isolate,
        records=refseq_records,
    )

    assert mock_isolate.accessions == {"MF062136"} | mock_isolate_accessions

    otu_after_promote = empty_repo.get_otu(otu_id)

    assert (
        otu_after_promote.accessions
        == new_refseq_isolate.accessions | mock_isolate_accessions - {"MF062136"}
    )

    assert otu_after_promote.get_isolate(isolate_init.id).accessions == (
        mock_isolate_accessions - {"MF062136"} | {"NC_055390"}
    )
