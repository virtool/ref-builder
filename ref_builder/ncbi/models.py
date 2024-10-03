"""Models for NCBI Genbank and Taxonomy data."""

from enum import StrEnum
from typing import Annotated, Any

from pydantic import (
    AliasChoices,
    BaseModel,
    ConfigDict,
    Field,
    computed_field,
    field_validator,
    model_validator,
)

from ref_builder.models import MolType, Strandedness, Topology


class NCBIDatabase(StrEnum):
    """NCBI databases used in ref-builder."""

    NUCCORE = "nuccore"
    TAXONOMY = "taxonomy"


class NCBIRank(StrEnum):
    """NCBI ranks used in ref-builder."""

    SPECIES = "species"
    ISOLATE = "isolate"


class NCBISourceMolType(StrEnum):
    """In vivo molecule types.

    Based on the INSDC controlled vocabulary list for the /mol_type qualifier

    Reference:
    https://www.insdc.org/submitting-standards/controlled-vocabulary-moltype-qualifier/
    """

    GENOMIC_DNA = "genomic DNA"
    OTHER_DNA = "other DNA"
    UNASSIGNED_DNA = "unassigned DNA"
    GENOMIC_RNA = "genomic RNA"
    MRNA = "mRNA"
    TRNA = "tRNA"
    TRANSCRIBED_RNA = "transcribed RNA"
    VIRAL_CRNA = "viral cRNA"
    OTHER_RNA = "other RNA"


class NCBISource(BaseModel):
    """An NCBI source table."""

    taxid: int
    organism: str
    mol_type: NCBISourceMolType
    isolate: str = ""
    host: str = ""
    segment: str = ""
    strain: str = ""
    clone: str = ""
    proviral: bool = False
    macronuclear: bool = False
    focus: bool = False
    transgenic: bool = False

    @model_validator(mode="before")
    @classmethod
    def db_xref_to_taxid(cls, data: dict) -> dict:
        """Parse db_xref if ``taxid`` is not provided to the model directly.

        This is the realistic use case. The source table does not contain a taxid field,
        but we want to pass a fake value to the model in testing.

        """
        if data.get("taxid"):
            return data

        if db_xref := data.get("db_xref"):
            data["taxid"] = db_xref.split(":")[1]
            return data

        raise ValueError("No db_xref or taxid value found in source table")


class NCBIGenbank(BaseModel):
    """An NCBI Genbank record."""

    model_config = ConfigDict(populate_by_name=True)

    accession: Annotated[str, Field(validation_alias="GBSeq_primary-accession")]
    accession_version: Annotated[str, Field(validation_alias="GBSeq_accession-version")]
    strandedness: Annotated[
        Strandedness,
        Field(validation_alias="GBSeq_strandedness"),
    ]
    moltype: Annotated[MolType, Field(validation_alias="GBSeq_moltype")]
    topology: Annotated[Topology, Field(validation_alias="GBSeq_topology")]
    definition: Annotated[str, Field(validation_alias="GBSeq_definition")]
    organism: Annotated[str, Field(validation_alias="GBSeq_organism")]
    sequence: Annotated[
        str,
        Field(
            validation_alias="GBSeq_sequence",
            pattern=r"^[ATCGRYKMSWBDHVNatcgrykmswbdhvn]+$",
        ),
    ]
    source: Annotated[NCBISource, Field(validation_alias="GBSeq_feature-table")]
    comment: Annotated[str, Field(validation_alias="GBSeq_comment")] = ""

    @computed_field()
    def refseq(self) -> bool:
        """Whether this is a RefSeq record.

        RefSeq records have accessions that start with "NC_".
        """
        return self.accession.startswith("NC_")

    @field_validator("sequence", mode="after")
    @classmethod
    def uppercase_sequence(cls, raw: str) -> str:
        """Force the sequence field to uppercase."""
        return raw.upper()

    @field_validator("source", mode="before")
    @classmethod
    def convert_source(cls, raw: NCBISource | list[dict[str:Any]]) -> NCBISource:
        """If the source field isn't a ``NCBISource`` object, convert it."""
        if isinstance(raw, NCBISource):
            return raw

        for feature in raw:
            if feature["GBFeature_key"] == "source":
                data = {}

                for qualifier in feature["GBFeature_quals"]:
                    data[qualifier["GBQualifier_name"]] = qualifier.get(
                        "GBQualifier_value",
                        True,
                    )

                return NCBISource(**data)

        raise ValueError("Feature table contains no ``source`` table.")

    @model_validator(mode="after")
    def check_source(self) -> "NCBIGenbank":
        """Check that the source organism matches the record organism."""
        if self.source.organism == self.organism:
            return self

        raise ValueError("Non-matching organism fields on record and source")


class NCBILineage(BaseModel):
    """An NCBI lineage record."""

    model_config = ConfigDict(populate_by_name=True)

    id: Annotated[int, Field(validation_alias="TaxId")]
    name: Annotated[str, Field(validation_alias="ScientificName")]
    rank: Annotated[str, Field(validation_alias="Rank")]


class NCBITaxonomyOtherNames(BaseModel):
    """An NCBI taxonomy record's other names."""

    acronym: Annotated[list[str], Field(validation_alias="Acronym")] = []
    genbank_acronym: Annotated[list[str], Field(validation_alias="GenbankAcronym")] = []
    equivalent_name: Annotated[list[str], Field(validation_alias="EquivalentName")] = []
    synonym: Annotated[list[str], Field(validation_alias="Synonym")] = []
    includes: Annotated[list[str], Field(validation_alias="Includes")] = []


class NCBITaxonomy(BaseModel):
    """An NCBI taxonomy record."""

    model_config = ConfigDict(populate_by_name=True)

    id: Annotated[int, Field(validation_alias="TaxId")]
    name: Annotated[str, Field(validation_alias="ScientificName")]
    other_names: Annotated[
        NCBITaxonomyOtherNames,
        Field(validation_alias="OtherNames"),
    ] = NCBITaxonomyOtherNames()
    lineage: Annotated[list[NCBILineage], Field(validation_alias="LineageEx")]
    rank: Annotated[NCBIRank, Field(validation_alias=AliasChoices("rank", "Rank"))]

    @computed_field
    def species(self) -> NCBILineage:
        """Return the species level taxon in the lineage."""
        if self.rank is NCBIRank.SPECIES:
            return NCBILineage(id=self.id, name=self.name, rank=self.rank)

        for item in self.lineage:
            if item.rank == "species":
                return item

        raise ValueError("No species level taxon found in lineage")
