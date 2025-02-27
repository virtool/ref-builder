from pathlib import Path
from uuid import UUID, uuid4

import orjson
import pytest

from ref_builder.models import Molecule, MolType, Strandedness, Topology
from ref_builder.plan import Plan, Segment, SegmentRule
from ref_builder.repo import Repo, GITIGNORE_CONTENTS
from ref_builder.resources import (
    RepoIsolate,
    RepoOTU,
    RepoSequence,
)
from ref_builder.utils import Accession, DataType, IsolateName, IsolateNameType


@pytest.fixture()
def initialized_repo(empty_repo: Repo):
    """Return a pre-initialized mock Repo."""
    with empty_repo.lock(), empty_repo.use_transaction():
        otu = empty_repo.create_otu(
            acronym="TMV",
            legacy_id=None,
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=Plan.new(
                [
                    Segment.new(
                        length=150,
                        length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                        name=None,
                    )
                ]
            ),
            taxid=12242,
        )

        sequence_1 = empty_repo.create_sequence(
            otu.id,
            "TMVABC.1",
            "TMV",
            None,
            otu.plan.segments[0].id,
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
    with empty_repo.lock(), empty_repo.use_transaction():
        return empty_repo.create_otu(
            acronym="TMV",
            legacy_id="abcd1234",
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=Plan.new(
                segments=[
                    Segment.new(
                        length=150,
                        length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                        name=None,
                    )
                ]
            ),
            taxid=12242,
        )


class TestNew:
    def test_ok(self, empty_repo: Repo, tmp_path: Path):
        """Test that creating a new ``Repo`` object returns the expected object and
        creates the expected directory structure.
        """
        assert empty_repo.path == tmp_path / "test_repo"
        assert empty_repo.last_id == 1

        assert empty_repo.meta.data_type == DataType.GENOME
        assert empty_repo.meta.name == "Generic Viruses"
        assert empty_repo.meta.organism == "virus"

        assert empty_repo.settings.default_segment_length_tolerance == 0.03

        assert (empty_repo.path / ".gitignore").exists()

        with open(empty_repo.path / ".gitignore", "r") as f:
            assert f.read() == "\n".join(GITIGNORE_CONTENTS) + "\n"

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
        plan = Plan.new(
            segments=[
                Segment.new(
                    length=150,
                    length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                    name=None,
                )
            ]
        )

        with empty_repo.lock(), empty_repo.use_transaction():
            otu = empty_repo.create_otu(
                acronym="TMV",
                legacy_id="abcd1234",
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                plan=plan,
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
            representative_isolate=None,
            plan=Plan(
                id=plan.id,
                segments=[
                    Segment(
                        id=plan.segments[0].id,
                        length=150,
                        length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                        name=None,
                        rule=SegmentRule.REQUIRED,
                    )
                ],
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
                    "id": str(plan.id),
                    "segments": [
                        {
                            "id": str(plan.segments[0].id),
                            "length": 150,
                            "length_tolerance": empty_repo.settings.default_segment_length_tolerance,
                            "name": None,
                            "rule": SegmentRule.REQUIRED,
                        }
                    ],
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
        with empty_repo.lock():
            with empty_repo.use_transaction():
                empty_repo.create_otu(
                    acronym="TMV",
                    legacy_id=None,
                    molecule=Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MolType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    name="Tobacco mosaic virus",
                    plan=Plan.new(
                        segments=[
                            Segment.new(
                                length=150,
                                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                                name=None,
                            )
                        ]
                    ),
                    taxid=12242,
                )

            with (
                empty_repo.use_transaction(),
                pytest.raises(
                    ValueError,
                    match="An OTU with the name 'Tobacco mosaic virus' already exists",
                ),
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
                    plan=Plan.new(
                        segments=[
                            Segment.new(
                                length=150,
                                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                                name=None,
                            )
                        ]
                    ),
                    taxid=438782,
                )

    def test_duplicate_legacy_id(self, empty_repo: Repo):
        """Test that creating an OTU with a legacy ID that already exists raises a
        ``ValueError``.
        """
        with empty_repo.lock(), empty_repo.use_transaction():
            empty_repo.create_otu(
                acronym="TMV",
                legacy_id="abcd1234",
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                plan=Plan.new(
                    segments=[
                        Segment.new(
                            length=150,
                            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                            name=None,
                        )
                    ]
                ),
                taxid=12242,
            )

        with (
            pytest.raises(
                ValueError,
                match="An OTU with the legacy ID 'abcd1234' already exists",
            ),
            empty_repo.lock(),
            empty_repo.use_transaction(),
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
                plan=Plan.new(
                    segments=[
                        Segment.new(
                            length=150,
                            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                            name=None,
                        )
                    ]
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

        with empty_repo.lock(), empty_repo.use_transaction():
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

        with empty_repo.lock(), empty_repo.use_transaction():
            empty_repo.create_isolate(
                otu.id,
                None,
                IsolateName(IsolateNameType.ISOLATE, "A"),
            )

        with (
            empty_repo.lock(),
            empty_repo.use_transaction(),
            pytest.raises(
                ValueError,
                match="Isolate name already exists: Isolate A",
            ),
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

        with empty_repo.lock(), empty_repo.use_transaction():
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

    with empty_repo.lock(), empty_repo.use_transaction():
        sequence = empty_repo.create_sequence(
            otu.id,
            "TMVABC.1",
            "TMV",
            None,
            otu.plan.segments[0].id,
            "ACGT",
        )

    assert sequence is not None

    assert sequence == RepoSequence(
        id=sequence.id,
        accession=Accession(key="TMVABC", version=1),
        definition="TMV",
        legacy_id=None,
        segment=otu.plan.segments[0].id,
        sequence="ACGT",
    )

    with open(empty_repo.path.joinpath("src", "00000003.json")) as f:
        event = orjson.loads(f.read())

    del event["timestamp"]

    assert event == {
        "data": {
            "id": str(sequence.id),
            "accession": {"key": "TMVABC", "version": 1},
            "definition": "TMV",
            "legacy_id": None,
            "segment": str(otu.plan.segments[0].id),
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
        """Test that getting an OTU returns the expected ``RepoOTU`` object including
        two isolates with one sequence each.
        """
        monopartite_plan = Plan.new(
            segments=[
                Segment.new(
                    length=150,
                    length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                    name=None,
                )
            ]
        )

        with empty_repo.lock():
            with empty_repo.use_transaction():
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

            segment_id = otu.plan.segments[0].id

            with empty_repo.use_transaction():
                sequence_1 = empty_repo.create_sequence(
                    otu.id,
                    "TMVABC.1",
                    "TMV",
                    None,
                    segment_id,
                    "ACGT",
                )

                isolate_a = empty_repo.create_isolate(
                    otu.id,
                    None,
                    IsolateName(IsolateNameType.ISOLATE, "A"),
                )

                empty_repo.link_sequence(
                    otu.id, isolate_a.id, sequence_id=sequence_1.id
                )

                sequence_2 = empty_repo.create_sequence(
                    otu.id,
                    "TMVABCB.1",
                    "TMV",
                    None,
                    segment_id,
                    "ACGTGGAGAGACC",
                )

                isolate_b = empty_repo.create_isolate(
                    otu.id,
                    None,
                    IsolateName(IsolateNameType.ISOLATE, "B"),
                )

                empty_repo.link_sequence(
                    otu.id, isolate_b.id, sequence_id=sequence_2.id
                )

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
                        segment=segment_id,
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
                        segment=segment_id,
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
                representative_isolate=None,
                plan=Plan(
                    id=monopartite_plan.id,
                    segments=[
                        Segment(
                            id=segment_id,
                            length=150,
                            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                            name=None,
                            rule=SegmentRule.REQUIRED,
                        )
                    ],
                ),
                taxid=12242,
                isolates=otu_contents,
            ).model_dump()
        )

        assert empty_repo.last_id == 8

    def test_retrieve_nonexistent_otu(self, initialized_repo: Repo):
        """Test that getting an OTU that does not exist returns ``None``."""
        assert initialized_repo.get_otu(uuid4()) is None

    def test_accessions(self, initialized_repo: Repo):
        """Test that the `accessions` property returns the expected accessions."""
        assert next(initialized_repo.iter_otus()).accessions == {"TMVABC"}

    def test_blocked_accessions(self, initialized_repo: Repo):
        """Test that the `blocked_accessions` property returns the expected set of
        accessions.
        """
        otu = initialized_repo.get_otu_by_taxid(12242)

        with initialized_repo.lock(), initialized_repo.use_transaction():
            sequence = initialized_repo.create_sequence(
                otu.id,
                "TMVABCB.1",
                "TMV",
                None,
                otu.plan.segments[0].id,
                "ACGTGGAGAGACC",
            )

            isolate = initialized_repo.create_isolate(
                otu.id,
                None,
                IsolateName(type=IsolateNameType.ISOLATE, value="B"),
            )

            initialized_repo.link_sequence(otu.id, isolate.id, sequence_id=sequence.id)

            initialized_repo.exclude_accession(otu.id, "GROK")
            initialized_repo.exclude_accession(otu.id, "TOK")

        assert initialized_repo.get_otu(otu.id).blocked_accessions == {
            "TMVABC",
            "TMVABCB",
            "GROK",
            "TOK",
        }


def test_get_otu_id_from_isolate_id(initialized_repo: Repo):
    """Test that the OTU id can be retrieved from a isolate ID contained within."""
    otu = next(initialized_repo.iter_otus())

    isolate = otu.isolates[0]

    assert initialized_repo.get_otu_id_by_isolate_id(isolate.id) == otu.id


class TestGetIsolate:
    def test_by_id(self, initialized_repo: Repo):
        """Test that getting an isolate returns the expected ``RepoIsolate`` object."""
        otu = next(initialized_repo.iter_otus())

        for isolate in otu.isolates:
            assert otu.get_isolate(isolate.id) in otu.isolates

    def test_by_name(self, initialized_repo: Repo):
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

        with initialized_repo.lock():
            with initialized_repo.use_transaction():
                sequence = initialized_repo.create_sequence(
                    otu.id,
                    accession="EMPTY1.1",
                    definition="TMV B",
                    legacy_id=None,
                    segment=otu.plan.segments[0].id,
                    sequence="GACCACGTGGAGA",
                )

                isolate_unnamed = initialized_repo.create_isolate(
                    otu.id,
                    None,
                    None,
                )

                initialized_repo.link_sequence(
                    otu.id, isolate_unnamed.id, sequence_id=sequence.id
                )

            otu_after = next(initialized_repo.iter_otus())

        assert len(otu_after.isolate_ids) == len(otu.isolate_ids) + 1
        assert otu_after.accessions == {"TMVABC", "EMPTY1"}
        assert otu_after.isolate_ids == {isolate_before.id, isolate_unnamed.id}

        isolate_unnamed_after = otu_after.get_isolate(isolate_unnamed.id)

        assert isolate_unnamed_after.name is None
        assert isolate_unnamed_after.accessions == {"EMPTY1"}


def test_get_isolate_id_from_partial(initialized_repo: Repo):
    """Test that an isolate id can be retrieved from a truncated ``partial`` string."""
    otu = next(initialized_repo.iter_otus())

    isolate = otu.isolates[0]

    assert initialized_repo.get_isolate_id_by_partial(str(isolate.id)[:8]) == isolate.id


class TestExcludeAccessions:
    """Test that accessions can be excluded from future fetches."""

    def test_ok(self, empty_repo: Repo):
        otu = init_otu(empty_repo)

        assert otu.excluded_accessions == set()

        accessions = {"TM100021", "TM100022", "TM100023"}

        with empty_repo.lock(), empty_repo.use_transaction():
            empty_repo.exclude_accessions(otu.id, accessions)

        assert empty_repo.get_otu(otu.id).excluded_accessions == accessions
        assert empty_repo.last_id == 3

        with open(empty_repo.path / "src" / f"{empty_repo.last_id:08}.json") as f:
            event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "accessions": ["TM100021", "TM100022", "TM100023"],
                "action": "exclude",
            },
            "id": 3,
            "query": {
                "otu_id": str(otu.id),
            },
            "type": "UpdateExcludedAccessions",
        }

        with empty_repo.lock(), empty_repo.use_transaction():
            empty_repo.exclude_accessions(otu.id, {"TM100024"})

        assert empty_repo.last_id == 4
        assert empty_repo.get_otu(otu.id).excluded_accessions == accessions | {
            "TM100024"
        }

    def test_exists(self, initialized_repo: Repo):
        """Test that an accession the exists in an OTU sequence cannot be added to the
        excluded accessions.
        """
        assert initialized_repo.last_id == 5

        with initialized_repo.lock():
            otu = initialized_repo.get_otu_by_taxid(12242)

            assert otu.excluded_accessions == set()

            accession = next(iter(otu.accessions))

            with initialized_repo.use_transaction():
                initialized_repo.exclude_accessions(otu.id, {f"{accession}.1"})

            assert initialized_repo.last_id == 5
            assert initialized_repo.get_otu_by_taxid(12242).excluded_accessions == set()

        with open(
            initialized_repo.path / "src" / f"{initialized_repo.last_id:08}.json"
        ) as f:
            event = orjson.loads(f.read())

        assert event["type"] != "UpdateExcludedAccessions"

    def test_already_excluded(self, empty_repo: Repo):
        """Test that an attempted redundant exclusion does not create a new `event."""
        otu_id = init_otu(empty_repo).id

        accessions = {"TM100021", "TM100022", "TM100023"}

        with empty_repo.lock():
            assert empty_repo.get_otu(otu_id).excluded_accessions == set()

            with empty_repo.use_transaction():
                empty_repo.exclude_accessions(otu_id, accessions)

            assert empty_repo.last_id == 3

            with empty_repo.use_transaction():
                empty_repo.exclude_accessions(otu_id, {"TM100021"})

            assert empty_repo.get_otu(otu_id).excluded_accessions == accessions

        assert empty_repo.last_id == 3

    def test_partially_already_excluded(self, empty_repo: Repo):
        """Test that a partially redundant list of exclusions creates a new event with
        the redundant accession omitted.
        """
        otu_id = init_otu(empty_repo).id

        accessions = {"TM100021", "TM100022", "TM100023"}

        with empty_repo.lock():
            with empty_repo.use_transaction():
                empty_repo.exclude_accessions(otu_id, accessions)

            otu_before = empty_repo.get_otu(otu_id)

            assert empty_repo.last_id == 3
            assert otu_before.excluded_accessions == accessions

            with empty_repo.use_transaction():
                empty_repo.exclude_accessions(otu_id, {"TM100023", "TM100024"})

            otu_after = empty_repo.get_otu(otu_id)

            assert empty_repo.last_id == 4
            assert otu_after.excluded_accessions == otu_before.excluded_accessions | {
                "TM100024"
            }

        with open(
            empty_repo.path.joinpath("src", f"{empty_repo.last_id:08}.json")
        ) as f:
            event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "accessions": ["TM100024"],
                "action": "exclude",
            },
            "id": 4,
            "query": {
                "otu_id": str(otu_id),
            },
            "type": "UpdateExcludedAccessions",
        }


class TestAllowAccessions:
    """Test that accessions allowed back into the OTU are no longer contained
    in the excluded accessions set.
    """

    def test_ok(self, empty_repo: Repo):
        """Test that Repo.allow_accessions() produces the correct event and
        creates the expected RepoOTU.excluded_accessions set.
        """
        otu = init_otu(empty_repo)

        accessions = {"TM100021", "TM100022", "TM100023"}

        with empty_repo.lock():
            with empty_repo.use_transaction():
                empty_repo.exclude_accessions(otu.id, accessions)

            assert empty_repo.get_otu(otu.id).excluded_accessions == accessions

            with empty_repo.use_transaction():
                empty_repo.allow_accessions(otu.id, ["TM100021", "TM100022"])

            otu_after = empty_repo.get_otu(otu.id)

        with open(empty_repo.path.joinpath("src", "00000004.json")) as f:
            event = orjson.loads(f.read())

            del event["timestamp"]

            assert event == {
                "data": {
                    "accessions": ["TM100021", "TM100022"],
                    "action": "allow",
                },
                "id": 4,
                "query": {
                    "otu_id": str(otu_after.id),
                },
                "type": "UpdateExcludedAccessions",
            }

        assert otu_after.excluded_accessions == {"TM100023"}

    def test_skip_redundant_accessions(self, empty_repo: Repo):
        """Test that an event only gets written if the accession exists
        in the exclusion list, avoiding the creation of a redundant event.
        """
        otu = init_otu(empty_repo)

        accessions = {"TM100021", "TM100022", "TM100023"}

        with empty_repo.lock():
            with empty_repo.use_transaction():
                empty_repo.exclude_accessions(otu.id, accessions)

            assert empty_repo.get_otu(otu.id).excluded_accessions == accessions

            with empty_repo.use_transaction():
                empty_repo.allow_accessions(otu.id, ["TM100024"])

            assert empty_repo.get_otu(otu.id).excluded_accessions == accessions

            with empty_repo.use_transaction():
                empty_repo.allow_accessions(otu.id, accessions)

            otu_after = empty_repo.get_otu(otu.id).excluded_accessions == set()

        assert not empty_repo.path.joinpath("src", "00000005.json").exists()


def test_delete_isolate(initialized_repo: Repo):
    """Test that an isolate can be redacted from an OTU."""
    with initialized_repo.lock():
        otu_before = next(initialized_repo.iter_otus())

        otu_id = otu_before.id

        isolate_before = next(iter(otu_before.isolates))

        with initialized_repo.use_transaction():
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


def test_replace_sequence(initialized_repo: Repo):
    """Test the replacement of am existing sequence using a new Genbank accession and
    record.
    """
    otu_before = initialized_repo.get_otu_by_taxid(12242)

    accession = "TMVABC"

    isolate_id, replaced_sequence_id = (
        otu_before.get_sequence_id_hierarchy_from_accession(accession)
    )

    with initialized_repo.lock():
        with initialized_repo.use_transaction():
            new_sequence = initialized_repo.create_sequence(
                otu_before.id,
                "TMVABCC.1",
                "TMV edit",
                None,
                otu_before.plan.segments[0].id,
                "ACGTGGAGAGACCA",
            )

        with initialized_repo.use_transaction():
            initialized_repo.replace_sequence(
                otu_id=otu_before.id,
                isolate_id=isolate_id,
                sequence_id=new_sequence.id,
                replaced_sequence_id=replaced_sequence_id,
                rationale="Testing sequence redaction",
            )

        otu_after = initialized_repo.get_otu(otu_before.id)

    assert otu_after.get_sequence_by_id(new_sequence.id) is not None
    assert otu_after.get_sequence_by_id(replaced_sequence_id) is None


class TestMalformedEvent:
    """Test that malformed events cannot be rehydrated."""

    def test_bad_event_typing(self, initialized_repo: Repo):
        """Test that an event with an invalid event type discriminator does not attempt
        to rehydrate.
        """
        filepath = initialized_repo.path.joinpath("src", "00000002.json")

        with open(filepath, "rb") as f:
            event = orjson.loads(f.read())

        otu = initialized_repo.get_otu_by_taxid(12242)

        assert type(otu) is RepoOTU

        event["type"] = "MalformedEvent"

        with open(filepath, "wb") as f:
            f.write(orjson.dumps(event))

        with pytest.raises(ValueError, match="Unknown event type: MalformedEvent"):
            initialized_repo.get_otu_by_taxid(12242)

    def test_bad_event_data(self, initialized_repo: Repo):
        """Test that an event with bad data cannot be rehydrated."""
        path = initialized_repo.path.joinpath("src", "00000002.json")

        with open(path, "rb") as f:
            event = orjson.loads(f.read())

        with initialized_repo.lock():
            otu = initialized_repo.get_otu_by_taxid(12242)

        assert type(otu) is RepoOTU

        event["data"]["taxid"] = "popcorn"

        with open(path, "wb") as f:
            f.write(orjson.dumps(event))

        with (
            initialized_repo.lock(),
            pytest.raises(
                ValueError,
                match="Input should be a valid integer, unable to parse string",
            ),
        ):
            initialized_repo.get_otu_by_taxid(12242)
