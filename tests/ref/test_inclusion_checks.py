import pytest
from pydantic import ValidationError
from faker import Faker

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.isolate import create_isolate
from ref_builder.otu.utils import IsolateName, IsolateNameType, check_sequence_length
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory
from tests.fixtures.providers import AccessionProvider, SequenceProvider

faker = Faker()
faker.add_provider(AccessionProvider)
faker.add_provider(SequenceProvider)


class TestSequenceLengthCheck:
    """Test sequence length check functions."""

    @pytest.mark.parametrize(
        ("sequence_lengths", "tolerance"),
        [
            ([1000, 970, 1030], 0.03),
            ([1000, 980, 1020], 0.02),
            ([1000, 960, 1040], 0.04),
        ],
    )
    def test_check_sequence_length_function_success(
        self, sequence_lengths: list[int], tolerance: float
    ):
        """Test that check_sequence_length() returns True
        given different sequences and tolerance levels.
        """
        for sequence_length in sequence_lengths:
            fake_sequence = faker.sequence(min=sequence_length, max=sequence_length)

            assert check_sequence_length(fake_sequence, 1000, tolerance)

    @pytest.mark.parametrize(
        ("sequence_lengths", "tolerance"),
        [
            ([970, 1030], 0.02),
            ([969, 1031], 0.03),
            ([979, 1021], 0.02),
            ([959, 1041], 0.04),
        ],
    )
    def test_check_sequence_length_function_fail(
        self, sequence_lengths: list[int], tolerance: float
    ):
        """Test that check_sequence_length() returns False
        given different sequences and tolerance levels.
        """
        for sequence_length in sequence_lengths:
            fake_sequence = faker.sequence(min=sequence_length, max=sequence_length)

            assert not check_sequence_length(fake_sequence, 1000, tolerance)


class TestAddMultipartiteIsolate:
    """Test that new multipartite isolates follow plan length limits."""

    @pytest.fixture(autouse=True)
    def _setup(
        self,
        ncbi_genbank_factory: NCBIGenbankFactory,
        ncbi_source_factory: NCBISourceFactory,
    ):
        def func(
            otu: RepoOTU,
            sequence_length_multiplier: float = 1.0,
        ) -> list[NCBIGenbank]:
            """Generate a collection of mock Genbank records capable of passing the isolate inclusion checks
            on a given multipartite OTU, as long as sequence_length_multiplier is 1.0.
            Builds one mock record per segment.

            :param otu: A monopartite OTU that the generated record must match.
            :param sequence_length_multiplier: A float multiplier for the generated mock sequence length.
            """
            records = []

            accession_starter = 100000

            for i in range(len(otu.plan.required_segments)):
                segment = otu.plan.required_segments[i]

                source = ncbi_source_factory.build(
                    taxid=otu.taxid,
                    organism=otu.name,
                    segment=str(segment.name),
                )

                sequence_length = int(segment.length * sequence_length_multiplier)

                records.append(
                    ncbi_genbank_factory.build(
                        accession=f"AB{accession_starter + i}",
                        sequence=faker.sequence(
                            min=sequence_length, max=sequence_length
                        ),
                        source=source,
                    )
                )

            return records

        self.create_mock_isolate_records = func

    @pytest.mark.parametrize("sequence_length_multiplier", [1.0, 1.03, 0.98])
    def test_ok(self, scratch_repo: Repo, sequence_length_multiplier: float):
        """Test that sequences within recommended length tolerance
        are added without issue.
        """
        otu_before = scratch_repo.get_otu_by_taxid(438782)

        records = self.create_mock_isolate_records(
            otu_before,
            sequence_length_multiplier,
        )

        with scratch_repo.lock(), scratch_repo.use_transaction():
            isolate = create_isolate(
                scratch_repo,
                otu_before,
                IsolateName(type=IsolateNameType.ISOLATE, value="mock"),
                records,
            )

        otu_after = scratch_repo.get_otu_by_taxid(438782)

        assert isolate.id in otu_after.isolate_ids
        assert isolate.accessions.issubset(otu_after.accessions)

    @pytest.mark.parametrize("sequence_length_multiplier", [0.5, 1.035, 20.0])
    def test_fail(self, scratch_repo: Repo, sequence_length_multiplier: float):
        """Test that sequences that exceed recommended length tolerance are
        automatically rejected.
        """
        otu_before = scratch_repo.get_otu_by_taxid(438782)

        records = self.create_mock_isolate_records(
            otu_before,
            sequence_length_multiplier,
        )

        with scratch_repo.lock(), scratch_repo.use_transaction():
            try:
                create_isolate(
                    scratch_repo,
                    otu_before,
                    IsolateName(type=IsolateNameType.ISOLATE, value="mock"),
                    records,
                )
            except ValidationError as exc:
                for error in exc.errors():
                    assert error["type"] in ("sequence_too_short", "sequence_too_long")

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert otu_after.isolate_ids == otu_before.isolate_ids
        assert otu_after.accessions == otu_before.accessions
