"""Tests for OTU models."""

import uuid

import pytest

from ref_builder.otu.models import OTU, OTUBase
from tests.fixtures.factories import OTUFactory


class TestOTU:
    """Test the ``OTU`` model which is used for complete validation of OTUs."""

    otu: OTUBase

    @pytest.fixture(autouse=True)
    def _build_otu(self, otu_factory: OTUFactory):
        self.otu = otu_factory.build()

    def test_ok(self):
        """Test that a valid OTU passes validation."""
        assert OTU.model_validate(self.otu.model_dump())

    def test_excluded_accessions(self):
        """Test that validation fails if the OTU includes accessions that are included
        in ``excluded_accessions``.
        """
        accession = self.otu.isolates[0].sequences[0].accession.key

        self.otu.excluded_accessions.add(accession)

        with pytest.raises(
            ValueError, match=f"Excluded accessions found in the OTU: {accession}"
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_no_isolates(self):
        """Test that validation fails if the OTU has no isolates."""
        self.otu.isolates = []
        self.otu.sequences = []

        with pytest.raises(
            ValueError, match="List should have at least 1 item after validation, not 0"
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_no_representative_isolate(self):
        """Test that validation fails if the OTU has no representative isolate."""
        self.otu.representative_isolate = None

        with pytest.raises(
            ValueError,
            match="UUID input should be a string, bytes or UUID object.*input_value=None",
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_missing_representative_isolate(self):
        """Test that validation fails if the representative isolate is not in the
        OTU.
        """
        self.otu.representative_isolate = uuid.uuid4()

        with pytest.raises(
            ValueError, match="Representative isolate must be in the OTU"
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_unique_isolate_names(self):
        """Test that validation fails if the OTU has duplicate isolate names."""
        duplicate_name = self.otu.isolates[1].name = self.otu.isolates[0].name

        with pytest.raises(
            ValueError,
            match=f"Isolate names must be unique. Non-unique names: {duplicate_name}",
        ):
            assert OTU.model_validate(self.otu.model_dump())
