import pytest
from syrupy import SnapshotAssertion
from syrupy.filters import props

from ref_builder.repo import Repo
from ref_builder.resources import (
    RepoIsolate,
    RepoOTU,
    RepoSequence,
)
from ref_builder.snapshotter.models import (
    OTUSnapshotIsolate,
    OTUSnapshotOTU,
    OTUSnapshotSequence,
)


def test_otu_model_adherence(scratch_repo: Repo):
    """Check the OTU snapshot model for missing fields relative to RepoOTU."""
    otu = scratch_repo.get_otu_by_taxid(1169032)

    otu_fields = set(otu.model_dump().keys())

    otu_fields.remove("isolates")
    otu_fields.remove("excluded_accessions")
    otu_fields.remove("schema")

    otu_fields.add("otu_schema")

    for field in otu_fields:
        assert field in OTUSnapshotOTU.model_fields


class TestRepoToSnapshotModel:
    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [
            (
                3158377,
                [
                    "NC_010314",
                    "NC_010315",
                    "NC_010316",
                    "NC_010317",
                    "NC_010318",
                    "NC_010319",
                ],
            ),
            (
                1169032,
                ["MH200607", "NC_003355", "KJ207375", "MK431779", "AB017504"],
            ),
        ],
    )
    def test_sequence_conversion(
        self,
        taxid: int,
        accessions: list[str],
        scratch_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for accession in accessions:
            original_sequence = otu.get_sequence_by_accession(accession)

            assert type(original_sequence) is RepoSequence

            converted_model = OTUSnapshotSequence(**original_sequence.dict())

            assert converted_model.model_dump() == snapshot(exclude=props("id"))

    @pytest.mark.parametrize("taxid", [1441799, 430059])
    def test_isolate_conversion(
        self,
        taxid: int,
        scratch_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        for isolate in otu.isolates:
            assert type(isolate) is RepoIsolate
            assert OTUSnapshotIsolate(**isolate.model_dump()).model_dump() == snapshot(
                exclude=props("id"),
            )

    @pytest.mark.parametrize("taxid", [1441799, 430059])
    def test_otu_conversion(
        self,
        taxid: int,
        scratch_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        otu = scratch_repo.get_otu_by_taxid(taxid)

        assert type(otu) is RepoOTU

        converted_model = OTUSnapshotOTU(**otu.model_dump())

        assert converted_model.id == otu.id

        assert converted_model.model_dump(by_alias=True) == snapshot(
            exclude=props("id", "repr_isolate"),
        )
