from pathlib import Path
from uuid import UUID, uuid4

import orjson
import pytest

from ref_builder.models import Molecule, MolType, Strandedness, Topology
from ref_builder.plan import MonopartitePlan
from ref_builder.repo import Repo
from ref_builder.resources import (
    RepoIsolate,
    RepoOTU,
    RepoSequence,
)
from ref_builder.utils import Accession, DataType, IsolateName, IsolateNameType


@pytest.fixture()
def initialized_repo(empty_repo: Repo):
    """Return a pre-initialized mock Repo."""
    otu = empty_repo.create_otu(
        acronym="TMV",
        legacy_id=None,
        molecule=Molecule(
            strandedness=Strandedness.SINGLE,
            type=MolType.RNA,
            topology=Topology.LINEAR,
        ),
        name="Tobacco mosaic virus",
        plan=MonopartitePlan.new(
            length=150,
            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
        ),
        taxid=12242,
    )

    sequence_1 = empty_repo.create_sequence(
        otu.id,
        "TMVABC.1",
        "TMV",
        None,
        "RNA",
        "ACGT",
    )

    isolate_a = empty_repo.create_isolate(
        otu.id,
        None,
        IsolateName(IsolateNameType.ISOLATE, "A"),
    )

    empty_repo.link_sequence(otu.id, isolate_a.id, sequence_1.id)

    return empty_repo


def init_otu(empty_repo: Repo) -> RepoOTU:
    return empty_repo.create_otu(
        acronym="TMV",
        legacy_id="abcd1234",
        molecule=Molecule(
            strandedness=Strandedness.SINGLE,
            type=MolType.RNA,
            topology=Topology.LINEAR,
        ),
        name="Tobacco mosaic virus",
        plan=MonopartitePlan.new(
            length=150,
            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
        ),
        taxid=12242,
    )


class TestNew:
    def test_ok(self, empty_repo: Repo, tmp_path: Path):
        """Test that creating a new ``Repo`` object returns the expected object and creates
        the expected directory structure.
        """
        assert empty_repo.path == tmp_path / "test_repo"
        assert empty_repo.last_id == 1

        assert empty_repo.meta.data_type == DataType.GENOME
        assert empty_repo.meta.name == "Generic Viruses"
        assert empty_repo.meta.organism == "virus"

        assert empty_repo.settings.default_segment_length_tolerance == 0.03

    def test_alternate_settings(self, tmp_path: Path):
        """Test retrieval of non-default settings."""
        repo = Repo.new(
            DataType.GENOME,
            "Generic Viruses",
            tmp_path / "alt_setting_repo",
            "virus",
            default_segment_length_tolerance=0.05,
        )

        assert repo.settings.default_segment_length_tolerance == 0.05


class TestCreateOTU:
    def test_ok(self, empty_repo: Repo):
        """Test that creating an OTU returns the expected ``RepoOTU`` object and creates
        the expected event file.
        """
        monopartite_plan = MonopartitePlan.new(
            length=150,
            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
        )

        otu = empty_repo.create_otu(
            acronym="TMV",
            legacy_id="abcd1234",
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=monopartite_plan,
            taxid=12242,
        )

        assert otu == RepoOTU(
            id=otu.id,
            acronym="TMV",
            excluded_accessions=set(),
            legacy_id="abcd1234",
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            repr_isolate=None,
            plan=MonopartitePlan(
                id=monopartite_plan.id,
                plan_type="monopartite",
                length=150,
                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
            ),
            taxid=12242,
            isolates=[],
        )

        with open(empty_repo.path.joinpath("src", "00000002.json")) as f:
            event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "id": str(otu.id),
                "acronym": "TMV",
                "molecule": {
                    "strandedness": "single",
                    "type": "RNA",
                    "topology": "linear",
                },
                "legacy_id": "abcd1234",
                "name": "Tobacco mosaic virus",
                "plan": {
                    "id": str(monopartite_plan.id),
                    "length": 150,
                    "length_tolerance": empty_repo.settings.default_segment_length_tolerance,
                    "name": None,
                    "plan_type": "monopartite",
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
        otu = empty_repo.create_otu(
            acronym="TMV",
            legacy_id=None,
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=MonopartitePlan.new(
                length=150,
                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
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
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                plan=MonopartitePlan.new(
                    length=150,
                    length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                ),
                taxid=438782,
            )

    def test_duplicate_legacy_id(self, empty_repo: Repo):
        """Test that creating an OTU with a legacy ID that already exists raises a
        ``ValueError``.
        """
        otu = empty_repo.create_otu(
            acronym="TMV",
            legacy_id="abcd1234",
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=MonopartitePlan.new(
                length=150,
                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
            ),
            taxid=12242,
        )

        with pytest.raises(
            ValueError,
            match="An OTU with the legacy ID 'abcd1234' already exists",
        ):
            empty_repo.create_otu(
                acronym="",
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                legacy_id="abcd1234",
                name="Abaca bunchy top virus",
                plan=MonopartitePlan.new(
                    length=150,
                    length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                ),
                taxid=438782,
            )


class TestCreateIsolate:
    """Test the creation and addition of new isolates in Repo."""

    def test_ok(self, empty_repo: Repo):
        """Test that creating an isolate returns the expected ``RepoIsolate`` object and
        creates the expected event file.
        """
        otu = init_otu(empty_repo)

        isolate = empty_repo.create_isolate(
            otu.id,
            None,
            IsolateName(IsolateNameType.ISOLATE, "A"),
        )

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

    def test_name_exists(self, empty_repo: Repo):
        """Test that a ValueError is raised if an isolate name is already taken."""
        otu = init_otu(empty_repo)

        empty_repo.create_isolate(
            otu.id,
            None,
            IsolateName(IsolateNameType.ISOLATE, "A"),
        )

        with pytest.raises(
            ValueError,
            match="Isolate name already exists: Isolate A",
        ):
            empty_repo.create_isolate(
                otu.id,
                None,
                IsolateName(IsolateNameType.ISOLATE, "A"),
            )

    def test_create_unnamed(self, empty_repo):
        """Test that creating an isolate returns the expected ``RepoIsolate`` object and
        creates the expected event file.
        """
        otu = init_otu(empty_repo)

        isolate = empty_repo.create_isolate(
            otu.id,
            None,
            None,
        )

        assert isinstance(isolate.id, UUID)
        assert isolate.sequences == []
        assert isolate.name is None


def test_create_sequence(empty_repo: Repo):
    """Test that creating a sequence returns the expected ``RepoSequence`` object and
    creates the expected event file.
    """
    otu = init_otu(empty_repo)

    sequence = empty_repo.create_sequence(
        otu.id,
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

    with open(empty_repo.path.joinpath("src", "00000003.json")) as f:
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
        "id": 3,
        "query": {
            "otu_id": str(otu.id),
            "sequence_id": str(sequence.id),
        },
        "type": "CreateSequence",
    }

    assert empty_repo.last_id == 3


class TestGetOTU:
    """Test the retrieval of OTU data."""

    def test_ok(self, empty_repo: Repo):
        """Test that getting an OTU returns the expected ``RepoOTU`` object including two
        isolates with one sequence each.
        """
        monopartite_plan = MonopartitePlan.new(
            length=150,
            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
        )

        otu = empty_repo.create_otu(
            acronym="TMV",
            legacy_id=None,
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            taxid=12242,
            plan=monopartite_plan,
        )

        sequence_1 = empty_repo.create_sequence(
            otu.id,
            "TMVABC.1",
            "TMV",
            None,
            "RNA",
            "ACGT",
        )

        isolate_a = empty_repo.create_isolate(
            otu.id,
            None,
            IsolateName(IsolateNameType.ISOLATE, "A"),
        )

        empty_repo.link_sequence(otu.id, isolate_a.id, sequence_id=sequence_1.id)

        sequence_2 = empty_repo.create_sequence(
            otu.id,
            "TMVABCB.1",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        isolate_b = empty_repo.create_isolate(
            otu.id,
            None,
            IsolateName(IsolateNameType.ISOLATE, "B"),
        )

        empty_repo.link_sequence(otu.id, isolate_b.id, sequence_id=sequence_2.id)

        otu = empty_repo.get_otu(otu.id)

        otu_contents = [
            RepoIsolate(
                id=isolate_a.id,
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
                id=isolate_b.id,
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
            otu.model_dump()
            == RepoOTU(
                id=otu.id,
                acronym="TMV",
                excluded_accessions=set(),
                legacy_id=None,
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                repr_isolate=None,
                plan=MonopartitePlan(
                    id=monopartite_plan.id,
                    plan_type="monopartite",
                    length=150,
                    length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                ),
                taxid=12242,
                isolates=otu_contents,
            ).model_dump()
        )

        assert empty_repo.last_id == 8

    def test_retrieve_nonexistent_otu(self, initialized_repo: Repo):
        """Test that getting an OTU that does not exist returns ``None``."""
        assert initialized_repo.get_otu(uuid4()) is None

    def test_get_accessions(self, initialized_repo: Repo):
        otu = next(initialized_repo.iter_otus())

        assert otu.accessions == {"TMVABC"}

        sequence_2 = initialized_repo.create_sequence(
            otu.id,
            "TMVABCB.1",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        isolate_b = initialized_repo.create_isolate(
            otu.id,
            None,
            IsolateName(type=IsolateNameType.ISOLATE, value="B"),
        )

        initialized_repo.link_sequence(otu.id, isolate_b.id, sequence_id=sequence_2.id)

        otu = next(initialized_repo.iter_otus())

        assert otu.accessions == {"TMVABC", "TMVABCB"}

    def test_get_blocked_accessions(self, initialized_repo: Repo):
        otu = initialized_repo.get_otu_by_taxid(12242)

        sequence_2 = initialized_repo.create_sequence(
            otu.id,
            "TMVABCB.1",
            "TMV",
            None,
            "RNA",
            "ACGTGGAGAGACC",
        )

        isolate_b = initialized_repo.create_isolate(
            otu.id,
            None,
            IsolateName(type=IsolateNameType.ISOLATE, value="B"),
        )

        initialized_repo.link_sequence(otu.id, isolate_b.id, sequence_id=sequence_2.id)

        initialized_repo.exclude_accession(otu.id, "GROK")
        initialized_repo.exclude_accession(otu.id, "TOK")

        otu = initialized_repo.get_otu(otu.id)

        assert otu.blocked_accessions == {"TMVABC", "TMVABCB", "GROK", "TOK"}


class TestGetIsolate:
    def test_get_isolate(self, initialized_repo: Repo):
        """Test that getting an isolate returns the expected ``RepoIsolate`` object."""
        otu = next(initialized_repo.iter_otus())

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

    def test_get_with_unnamed_isolate(self, initialized_repo: Repo):
        """Test that getting an OTU with an unnamed isolate ID behaves as expected."""
        otu = next(initialized_repo.iter_otus())

        isolate_before = otu.isolates[0]

        sequence_1_1 = initialized_repo.create_sequence(
            otu.id,
            accession="EMPTY1.1",
            definition="TMV B",
            legacy_id=None,
            segment="RNA A",
            sequence="GACCACGTGGAGA",
        )

        sequence_2_1 = initialized_repo.create_sequence(
            otu.id,
            accession="EMPTY2.1",
            definition="TMV A",
            legacy_id=None,
            segment="RNA B",
            sequence="ACTAAGAGAAAAA",
        )

        isolate_unnamed = initialized_repo.create_isolate(
            otu.id,
            None,
            None,
        )

        initialized_repo.link_sequence(
            otu.id, isolate_unnamed.id, sequence_id=sequence_1_1.id
        )

        initialized_repo.link_sequence(
            otu.id, isolate_unnamed.id, sequence_id=sequence_2_1.id
        )

        otu_after = next(initialized_repo.iter_otus())

        assert len(otu_after.isolate_ids) == len(otu.isolate_ids) + 1

        assert "EMPTY1" in otu_after.accessions

        assert "EMPTY2" in otu_after.accessions

        assert otu_after.isolate_ids == {isolate_before.id, isolate_unnamed.id}

        isolate_unnamed_after = otu_after.get_isolate(isolate_unnamed.id)

        assert isolate_unnamed_after.name is None

        assert isolate_unnamed_after.accessions == {"EMPTY1", "EMPTY2"}


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


def test_delete_isolate(initialized_repo):
    """Test that an isolate can be redacted from an OTU."""
    otu_before = next(initialized_repo.iter_otus())

    otu_id = otu_before.id

    isolate_before = list(otu_before.isolates)[0]

    initialized_repo.delete_isolate(
        otu_id,
        isolate_before.id,
        rationale="Testing redaction",
    )

    otu_after = initialized_repo.get_otu(otu_id)

    assert otu_before != otu_after

    assert len(otu_after.isolates) == len(otu_before.isolates) - 1

    assert isolate_before.id not in otu_after.isolate_ids

    assert isolate_before.accessions not in otu_after.accessions


def test_replace_sequence(initialized_repo):
    """Test the replacement of am existing sequence using a new Genbank accession and record."""
    otu_before = initialized_repo.get_otu_by_taxid(12242)

    accession = "TMVABC"

    isolate_id, replaced_sequence_id = (
        otu_before.get_sequence_id_hierarchy_from_accession(accession)
    )

    new_sequence = initialized_repo.create_sequence(
        otu_before.id,
        "TMVABCC.1",
        "TMV edit",
        None,
        "RNA",
        "ACGTGGAGAGACCA",
    )

    initialized_repo.replace_sequence(
        otu_id=otu_before.id,
        isolate_id=isolate_id,
        sequence_id=new_sequence.id,
        replaced_sequence_id=replaced_sequence_id,
        rationale="Testing sequence redaction",
    )

    otu_after = initialized_repo.get_otu(otu_before.id)

    assert new_sequence.id in otu_after.sequence_ids

    assert replaced_sequence_id not in otu_after.sequence_ids


class TestMalformedEvent:
    """Test that malformed events cannot be rehydrated."""

    def test_bad_event_typing(self, initialized_repo):
        """Test that an event with an invalid event type discriminator does not attempt to rehydrate"""
        filepath = initialized_repo.path.joinpath("src", "00000002.json")

        with open(filepath, "rb") as f:
            event = orjson.loads(f.read())

        otu = initialized_repo.get_otu_by_taxid(12242)

        assert type(otu) is RepoOTU

        event["type"] = "MalformedEvent"

        with open(filepath, "wb") as f:
            f.write(orjson.dumps(event))

        with pytest.raises(ValueError):
            initialized_repo.get_otu_by_taxid(12242)

    def test_bad_event_data(self, initialized_repo):
        """Test that an event with bad data cannot be rehydrated"""
        filepath = initialized_repo.path.joinpath("src", "00000002.json")

        with open(filepath, "rb") as f:
            event = orjson.loads(f.read())

        otu = initialized_repo.get_otu_by_taxid(12242)

        assert type(otu) is RepoOTU

        event["data"]["taxid"] = "popcorn"

        with open(filepath, "wb") as f:
            f.write(orjson.dumps(event))

        with pytest.raises(ValueError):
            initialized_repo.get_otu_by_taxid(12242)
