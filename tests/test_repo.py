import uuid
from pathlib import Path
from uuid import UUID

import orjson
import pytest

from ref_builder.models import Molecule, MolType, Strandedness, Topology
from ref_builder.repo import Repo
from ref_builder.resources import (
    RepoIsolate,
    RepoOTU,
    RepoSequence,
)
from ref_builder.schema import OTUSchema, Segment
from ref_builder.utils import Accession, DataType, IsolateName, IsolateNameType


@pytest.fixture()
def initialized_repo(empty_repo: Repo):
    otu = empty_repo.create_otu(
        "TMV",
        None,
        "Tobacco mosaic virus",
        OTUSchema(
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            segments=[Segment(name="A", length=150, required=True)],
        ),
        12242,
    )

    isolate_a = empty_repo.create_isolate(otu.id, None, "A", IsolateNameType.ISOLATE)
    empty_repo.create_sequence(
        otu.id,
        isolate_a.id,
        "TMVABC.1",
        "TMV",
        None,
        "RNA",
        "ACGT",
    )

    return empty_repo


def init_otu(empty_repo: Repo) -> RepoOTU:
    return empty_repo.create_otu(
        acronym="TMV",
        legacy_id="abcd1234",
        name="Tobacco mosaic virus",
        schema=OTUSchema(
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            segments=[Segment(name="A", required=True, length=100)],
        ),
        taxid=12242,
    )


def test_new(empty_repo: Repo, tmp_path: Path):
    """Test that creating a new ``Repo`` object returns the expected object and creates
    the expected directory structure.
    """
    assert empty_repo.path == tmp_path / "test_repo"
    assert empty_repo.last_id == 1

    assert empty_repo.meta.data_type == DataType.GENOME
    assert empty_repo.meta.name == "Generic Viruses"
    assert empty_repo.meta.organism == "virus"


class TestCreateOTU:
    def test_ok(self, empty_repo: Repo):
        """Test that creating an OTU returns the expected ``RepoOTU`` object and creates
        the expected event file.
        """
        otu = empty_repo.create_otu(
            acronym="TMV",
            legacy_id="abcd1234",
            name="Tobacco mosaic virus",
            schema=OTUSchema(
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                segments=[Segment(name="A", required=True, length=100)],
            ),
            taxid=12242,
        )

        assert (
            otu.dict()
            == RepoOTU(
                uuid=otu.id,
                acronym="TMV",
                excluded_accessions=None,
                legacy_id="abcd1234",
                name="Tobacco mosaic virus",
                schema=OTUSchema(
                    molecule=Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MolType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    segments=[Segment(name="A", required=True, length=100)],
                ),
                taxid=12242,
                isolates=[],
            ).dict()
        )

        with open(empty_repo.path.joinpath("src", "00000002.json")) as f:
            event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "id": str(otu.id),
                "acronym": "TMV",
                "legacy_id": "abcd1234",
                "name": "Tobacco mosaic virus",
                "schema": {
                    "molecule": {
                        "strandedness": "single",
                        "type": "RNA",
                        "topology": "linear",
                    },
                    "segments": [{"length": 100, "name": "A", "required": True}],
                    "multipartite": False,
                },
                "taxid": 12242,
            },
            "id": 2,
            "query": {
                "otu_id": str(otu.id),
            },
            "type": "CreateOTU",
        }

        assert empty_repo.last_id == 2

    def test_duplicate_name(self, empty_repo: Repo):
        """Test that creating an OTU with a name that already exists raises a
        ``ValueError``.
        """
        empty_repo.create_otu(
            acronym="TMV",
            legacy_id=None,
            name="Tobacco mosaic virus",
            schema=OTUSchema(
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                segments=[Segment(name="A", required=True)],
            ),
            taxid=12242,
        )

        with pytest.raises(
            ValueError,
            match="An OTU with the name 'Tobacco mosaic virus' already exists",
        ):
            empty_repo.create_otu(
                acronym="TMV",
                legacy_id=None,
                name="Tobacco mosaic virus",
                schema=OTUSchema(
                    molecule=Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MolType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    segments=[Segment(name="A", required=True)],
                ),
                taxid=438782,
            )

    def test_duplicate_legacy_id(self, empty_repo: Repo):
        """Test that creating an OTU with a legacy ID that already exists raises a
        ``ValueError``.
        """
        empty_repo.create_otu(
            acronym="TMV",
            legacy_id="abcd1234",
            name="Tobacco mosaic virus",
            schema=OTUSchema(
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                segments=[Segment(name="A", required=True)],
            ),
            taxid=12242,
        )

        with pytest.raises(
            ValueError,
            match="An OTU with the legacy ID 'abcd1234' already exists",
        ):
            empty_repo.create_otu(
                acronym="",
                legacy_id="abcd1234",
                name="Abaca bunchy top virus",
                schema=OTUSchema(
                    molecule=Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MolType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    segments=[Segment(name="A", required=True)],
                ),
                taxid=438782,
            )


def test_create_isolate(empty_repo: Repo):
    """Test that creating an isolate returns the expected ``RepoIsolate`` object and
    creates the expected event file.
    """
    otu = init_otu(empty_repo)

    isolate = empty_repo.create_isolate(otu.id, None, "A", IsolateNameType.ISOLATE)

    assert isinstance(isolate.id, UUID)
    assert isolate.sequences == []
    assert isolate.name.value == "A"
    assert isolate.name.type == "isolate"

    with open(empty_repo.path.joinpath("src", "00000003.json")) as f:
        event = orjson.loads(f.read())

    del event["timestamp"]

    assert event == {
        "data": {
            "id": str(isolate.id),
            "legacy_id": None,
            "name": {"type": "isolate", "value": "A"},
        },
        "id": 3,
        "query": {
            "otu_id": str(otu.id),
            "isolate_id": str(isolate.id),
        },
        "type": "CreateIsolate",
    }

    assert empty_repo.last_id == 3


def test_create_sequence(empty_repo: Repo):
    """Test that creating a sequence returns the expected ``RepoSequence`` object and
    creates the expected event file.
    """
    otu = init_otu(empty_repo)

    isolate = empty_repo.create_isolate(otu.id, None, "A", IsolateNameType.ISOLATE)

    sequence = empty_repo.create_sequence(
        otu.id,
        isolate.id,
        "TMVABC.1",
        "TMV",
        None,
        "RNA",
        "ACGT",
    )

    assert sequence is not None

    assert sequence == RepoSequence(
        id=sequence.id,
        accession=Accession(key="TMVABC", version=1),
        definition="TMV",
        legacy_id=None,
        segment="RNA",
        sequence="ACGT",
    )

    with open(empty_repo.path.joinpath("src", "00000004.json")) as f:
        event = orjson.loads(f.read())

    del event["timestamp"]

    assert event == {
        "data": {
            "id": str(sequence.id),
            "accession": ["TMVABC", 1],
            "definition": "TMV",
            "legacy_id": None,
            "segment": "RNA",
            "sequence": "ACGT",
        },
        "id": 4,
        "query": {
            "otu_id": str(otu.id),
            "isolate_id": str(isolate.id),
            "sequence_id": str(sequence.id),
        },
        "type": "CreateSequence",
    }

    assert empty_repo.last_id == 4


class TestRetrieveOTU:
    def test_get_otu(self, empty_repo: Repo):
        """Test that getting an OTU returns the expected ``RepoOTU`` object including two
        isolates with one sequence each.
        """
        otu = empty_repo.create_otu(
            acronym="TMV",
            legacy_id=None,
            name="Tobacco mosaic virus",
            taxid=12242,
            schema=OTUSchema(
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                segments=[Segment(name="A", required=True)],
            ),
        )

        isolate_a = empty_repo.create_isolate(
            otu.id,
            None,
            "A",
            IsolateNameType.ISOLATE,
        )
        empty_repo.create_sequence(
            otu.id,
            isolate_a.id,
            "TMVABC.1",
            "TMV",
            None,
            "RNA",
            "ACGT",
        )

        isolate_b = empty_repo.create_isolate(
            otu.id,
            None,
            "B",
            IsolateNameType.ISOLATE,
        )
        empty_repo.create_sequence(
            otu.id,
            isolate_b.id,
            "TMVABCB.1",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        otu = empty_repo.get_otu(otu.id)

        otu_contents = [
            RepoIsolate(
                uuid=isolate_a.id,
                legacy_id=None,
                name=IsolateName(type=IsolateNameType.ISOLATE, value="A"),
                sequences=[
                    RepoSequence(
                        id=otu.isolates[0].sequences[0].id,
                        accession=Accession(key="TMVABC", version=1),
                        definition="TMV",
                        legacy_id=None,
                        segment="RNA",
                        sequence="ACGT",
                    ),
                ],
            ),
            RepoIsolate(
                uuid=isolate_b.id,
                legacy_id=None,
                name=IsolateName(type=IsolateNameType.ISOLATE, value="B"),
                sequences=[
                    RepoSequence(
                        id=otu.isolates[1].sequences[0].id,
                        accession=Accession(key="TMVABCB", version=1),
                        definition="TMV",
                        legacy_id=None,
                        segment="RNA",
                        sequence="ACGTGGAGAGACC",
                    ),
                ],
            ),
        ]

        assert (
            otu.dict()
            == RepoOTU(
                uuid=otu.id,
                acronym="TMV",
                excluded_accessions=[],
                legacy_id=None,
                name="Tobacco mosaic virus",
                schema=OTUSchema(
                    molecule=Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MolType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    segments=[Segment(name="A", required=True)],
                ),
                taxid=12242,
                isolates=otu_contents,
            ).dict()
        )

        assert empty_repo.last_id == 6

    def test_retrieve_nonexistent_otu(self, initialized_repo: Repo):
        """Test that getting an OTU that does not exist returns ``None``."""
        assert initialized_repo.get_otu(uuid.uuid4()) is None

    def test_get_accessions(self, initialized_repo: Repo):
        otu = next(initialized_repo.iter_otus(ignore_cache=True))

        assert otu.accessions == {"TMVABC"}

        isolate_b = initialized_repo.create_isolate(
            otu.id,
            None,
            "B",
            IsolateNameType.ISOLATE,
        )
        initialized_repo.create_sequence(
            otu.id,
            isolate_b.id,
            "TMVABCB.1",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        otu = next(initialized_repo.iter_otus(ignore_cache=True))

        assert otu.accessions == {"TMVABC", "TMVABCB"}

    def test_get_blocked_accessions(self, initialized_repo: Repo):
        otu = initialized_repo.get_otu_by_taxid(12242)

        isolate_b = initialized_repo.create_isolate(
            otu.id,
            None,
            "B",
            IsolateNameType.ISOLATE,
        )

        initialized_repo.create_sequence(
            otu.id,
            isolate_b.id,
            "TMVABCB.1",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        initialized_repo.exclude_accession(otu.id, "GROK")
        initialized_repo.exclude_accession(otu.id, "TOK")

        otu = initialized_repo.get_otu(otu.id)

        assert otu.blocked_accessions == {"TMVABC", "TMVABCB", "GROK", "TOK"}


class TestGetIsolate:
    def test_get_isolate(self, initialized_repo: Repo):
        """Test that getting an isolate returns the expected ``RepoIsolate`` object."""
        otu = next(initialized_repo.iter_otus(ignore_cache=True))

        isolate_ids = {isolate.id for isolate in otu.isolates}

        for isolate_id in isolate_ids:
            assert otu.get_isolate(isolate_id) in otu.isolates

    def test_get_isolate_id_by_name(self, initialized_repo: Repo):
        """Test that getting an isolate ID by name returns the expected ID."""
        otu = next(initialized_repo.iter_otus())

        isolate_ids = {isolate.id for isolate in otu.isolates}

        assert (
            otu.get_isolate_id_by_name(
                IsolateName(type=IsolateNameType.ISOLATE, value="A"),
            )
            in isolate_ids
        )


def test_exclude_accession(empty_repo: Repo):
    """Test that excluding an accession from an OTU writes the expected event file and
    returns the expected OTU objects.
    """
    otu = init_otu(empty_repo)

    empty_repo.exclude_accession(otu.id, "TMVABC.1")

    with open(empty_repo.path.joinpath("src", "00000003.json")) as f:
        event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "accession": "TMVABC.1",
            },
            "id": 3,
            "query": {
                "otu_id": str(otu.id),
            },
            "type": "ExcludeAccession",
        }

    assert empty_repo.get_otu(otu.id).excluded_accessions == {
        "TMVABC.1",
    }

    empty_repo.exclude_accession(otu.id, "ABTV")

    assert empty_repo.get_otu(otu.id).excluded_accessions == {
        "TMVABC.1",
        "ABTV",
    }


class TestDirectDelete:
    def test_delete_isolate(self, initialized_repo: Repo):
        """Test that an isolate can be redacted from an OTU."""
        otu_before = next(initialized_repo.iter_otus())

        otu_id = otu_before.id

        isolate_before = list(otu_before.isolates)[0]

        initialized_repo.delete_isolate(
            otu_id, isolate_before.id, rationale="Testing redaction"
        )

        otu_after = initialized_repo.get_otu(otu_id)

        assert otu_before != otu_after

        assert len(otu_after.isolates) == len(otu_before.isolates) - 1

        assert isolate_before.id not in otu_after.isolate_ids

        assert isolate_before.accessions not in otu_after.accessions

    def test_replace_sequence(self, initialized_repo: Repo):
        """Test that a sequence can be redacted from an OTU."""
        otu_before = initialized_repo.get_otu_by_taxid(12242)

        otu_id = otu_before.id

        accession = "TMVABC"

        isolate_id, replaced_sequence_id = otu_before.get_sequence_id_hierarchy_from_accession(accession)

        assert otu_before.get_isolate(isolate_id).accessions == {"TMVABC"}

        new_sequence = initialized_repo.replace_sequence(
            otu_id,
            isolate_id,
            "TMVABCC.1",
            "TMV edit",
            None,
            "RNA",
            "ACGTGGAGAGACCA",
            replaced_sequence_id=replaced_sequence_id,
            rationale="Testing redaction",
        )

        otu_after = initialized_repo.get_otu(otu_id)

        assert otu_before != otu_after

        assert len(otu_after.accessions) == len(otu_before.accessions)

        assert replaced_sequence_id not in otu_after.sequence_ids

        assert new_sequence.id in otu_after.sequence_ids

        assert otu_after.get_isolate(isolate_id).accessions == {"TMVABCC"}

