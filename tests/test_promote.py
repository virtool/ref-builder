import datetime

import arrow
import pytest
from structlog.testing import capture_logs

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.create import create_otu_with_taxid
from ref_builder.otu.promote import (
    replace_otu_sequence_from_record,
    promote_otu_accessions_from_records,
    correct_sequences_in_otu,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoIsolate
from ref_builder.utils import Accession
from tests.fixtures.factories import IsolateFactory


def test_replace_sequence(empty_repo):
    """Test OTU sequence replacement."""

    with empty_repo.lock():
        otu_init = create_otu_with_taxid(
            empty_repo, 2164102, ["MF062125", "MF062126", "MF062127"], acronym=""
        )

    assert otu_init.accessions == {"MF062125", "MF062126", "MF062127"}

    initial_rep_isolate = otu_init.get_isolate(otu_init.representative_isolate)

    initial_rep_sequence_ids = [
        sequence.id for sequence in initial_rep_isolate.sequences
    ]

    refseq_records = NCBIClient(ignore_cache=False).fetch_genbank_records(
        ["NC_055390", "NC_055391", "NC_055392"]
    )

    record_by_replaceable_sequence_id = {
        initial_rep_sequence_ids[0]: refseq_records[0],
        initial_rep_sequence_ids[1]: refseq_records[1],
        initial_rep_sequence_ids[2]: refseq_records[2],
    }

    with empty_repo.lock():
        for sequence_id in record_by_replaceable_sequence_id:
            assert replace_otu_sequence_from_record(
                empty_repo,
                otu_init,
                sequence_id,
                replacement_record=record_by_replaceable_sequence_id[sequence_id],
                exclude_accession=True,
            )

    assert empty_repo.get_otu(otu_init.id).accessions == {
        "NC_055390",
        "NC_055391",
        "NC_055392",
    }


def test_multi_linked_promotion(empty_repo: Repo):
    """Test the promotion of a sequence that is linked to more than one isolate."""
    with empty_repo.lock():
        otu_init = create_otu_with_taxid(
            empty_repo, 2164102, ["MF062125", "MF062126", "MF062127"], acronym=""
        )

    assert otu_init.accessions == {"MF062125", "MF062126", "MF062127"}

    segment_l = otu_init.plan.get_segment_by_name_key("L")

    segment_l_rep_sequence = otu_init.get_sequence_by_accession("MF062125")

    otu_id = otu_init.id

    mock_isolate = RepoIsolate.model_validate(
        IsolateFactory.build_on_plan(otu_init.plan).model_dump()
    )

    for iterator in range(len(mock_isolate.sequences)):
        mock_sequence = mock_isolate.sequences[iterator]

        if mock_sequence.segment == segment_l.id:
            mock_isolate.sequences[iterator] = segment_l_rep_sequence.model_copy()

    with empty_repo.lock(), empty_repo.use_transaction():
        isolate_init = empty_repo.create_isolate(
            otu_id, legacy_id=None, name=mock_isolate.name
        )

        empty_repo.link_sequence(
            otu_id,
            isolate_init.id,
            segment_l_rep_sequence.id,
        )

        accession_counter = 1
        for mock_sequence in mock_isolate.sequences:
            # skip segment L
            if mock_sequence.segment == segment_l.id:
                continue

            sequence_init = empty_repo.create_sequence(
                otu_id,
                accession=f"FA00000{accession_counter}.1",
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

            accession_counter += 1

    otu_after_mock_isolate = empty_repo.get_otu(otu_id)

    assert isolate_init.id in otu_after_mock_isolate.isolate_ids

    assert otu_after_mock_isolate.accessions == {
        "MF062125",
        "MF062126",
        "MF062127",
        "FA000001",
        "FA000002",
    }

    refseq_records = NCBIClient(ignore_cache=False).fetch_genbank_records(
        ["NC_055390", "NC_055391", "NC_055392"],
    )

    with empty_repo.lock():
        assert promote_otu_accessions_from_records(
            empty_repo, otu_after_mock_isolate, refseq_records
        ) == {"NC_055390", "NC_055391", "NC_055392"}

    otu_after_promote = empty_repo.get_otu(otu_id)

    assert otu_after_promote.accessions == {
        "NC_055390",
        "NC_055391",
        "NC_055392",
        "FA000001",
        "FA000002",
    }

    assert otu_after_promote.get_isolate(
        otu_init.representative_isolate
    ).accessions == {"NC_055390", "NC_055391", "NC_055392"}

    assert otu_after_promote.get_isolate(isolate_init.id).accessions == (
        {"NC_055390", "FA000001", "FA000002"}
    )


@pytest.mark.ncbi()
class TestCorrectSequencesInOTU:
    """Test OTU-wide outdated sequence correction."""

    def test_ok(self, precached_repo: Repo):
        """Test a simple fetch and replace correction."""
        with precached_repo.lock():
            otu_init = create_otu_with_taxid(
                precached_repo,
                196375,
                ["NC_004452.1"],
                acronym="",
            )

        assert "NC_004452" in otu_init.accessions

        outdated_sequence = otu_init.get_sequence_by_accession("NC_004452")

        containing_isolate_id = next(iter(otu_init.get_isolate_ids_containing_sequence_id(outdated_sequence.id)))

        containing_isolate_init = otu_init.get_isolate(containing_isolate_id)

        assert outdated_sequence.accession.version == 1

        with precached_repo.lock():
            updated_sequence_id = next(iter(correct_sequences_in_otu(precached_repo, otu_init)))

        otu_after = precached_repo.get_otu(otu_init.id)

        assert containing_isolate_init.id in otu_after.isolate_ids

        containing_isolate_after = otu_after.get_isolate(containing_isolate_id)

        assert containing_isolate_after.sequence_ids != containing_isolate_init.sequence_ids

        assert containing_isolate_after.sequence_ids == {updated_sequence_id}

        updated_sequence = otu_after.get_sequence_by_id(updated_sequence_id)

        assert updated_sequence.accession.key == outdated_sequence.accession.key

        assert updated_sequence.accession.version > outdated_sequence.accession.version

    def test_with_future_date_limit(self, precached_repo: Repo):
        """Test that setting modification_date_start to a future date does returns no new sequences."""
        with precached_repo.lock():
            otu_init = create_otu_with_taxid(
                precached_repo,
                196375,
                ["NC_004452.1"],
                acronym="",
            )

        assert "NC_004452" in otu_init.accessions

        assert otu_init.get_sequence_by_accession("NC_004452").accession == (
            Accession(key="NC_004452", version=1)
        )

        with precached_repo.lock(), capture_logs() as captured_logs:
            correct_sequences_in_otu(
                precached_repo,
                otu_init,
                modification_date_start=arrow.utcnow().naive + datetime.timedelta(days=1),
            )

        assert any(log.get("event") == "All sequences are up to date." for log in captured_logs)
