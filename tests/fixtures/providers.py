"""Custom data providers for generating dummy data."""

from collections import OrderedDict

from faker.providers import BaseProvider

from tests.fixtures.organism import (
    ORGANISM_DESCRIPTOR_ADJECTIVES,
    ORGANISM_DESCRIPTOR_NOUNS,
    ORGANISM_HOSTS,
    ORGANISM_PART_AND_DESCRIPTORS,
    ORGANISM_PARTS,
    ORGANISM_TYPES,
    ORGANISM_VIRUSES,
)

UNCOMMON_NUCLEOTIDE_PROBABILITY = 0.01
"""Probability of an uncommon nucleotide appearing in a sequence."""

NUCLEOTIDE_PROBABILITIES = OrderedDict(
    [
        ("A", 0.245),
        ("T", 0.245),
        ("C", 0.245),
        ("G", 0.245),
        ("R", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("Y", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("K", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("M", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("S", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("W", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("B", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("D", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("H", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("V", UNCOMMON_NUCLEOTIDE_PROBABILITY),
        ("N", UNCOMMON_NUCLEOTIDE_PROBABILITY),
    ],
)
"""Probabilities of each nucleotide appearing in a sequence."""

MIN_SEQUENCE_LENGTH = 100
"""Minimum length of a sequence the should be generated."""

MAX_SEQUENCE_LENGTH = 1500
"""Maximum length of a sequence the should be generated."""


class AccessionProvider(BaseProvider):
    """Raw accession provider based on GenBank's guidelines for accession numbers."""

    def accession(self) -> str:
        """Return a pseudorandom accession number."""
        if self.random_int(0, 10) > 7:
            return self.genbank_accession()

        return self.refseq_accession()

    def accessions(self, count: int) -> list[str]:
        """Return a list of pseudorandom, consecutive accession numbers."""
        if self.random_int(0, 10) > 7:
            return self.genbank_accessions(count)

        return self.refseq_accessions(count)

    def genbank_accession(self) -> str:
        """Return a pseudorandom non-RefSeq accession number."""
        if self.random_int(0, 10) > 6:
            return self.bothify("?#####").upper()

        return self.bothify("??######").upper()

    def genbank_accessions(self, count: int) -> list[str]:
        """Return a list of pseudorandom, consecutive accession numbers."""
        # accessions = []

        if self.random_int(0, 10) > 6:
            prefix = self.random_uppercase_letter()

            ceiling = 99999

            digit_format = "{0:05d}"

        else:
            prefix = self.random_uppercase_letter() + self.random_uppercase_letter()

            ceiling = 999999

            digit_format = "{0:06d}"

        first_number = self.random_int(0, ceiling - count)

        return [
            f"{prefix}{digit_format.format(first_number + i)}" for i in range(count)
        ]

    def refseq_accession(self) -> str:
        """Return a pseudorandom RefSeq accession number."""
        return self.numerify("NC_######")

    def refseq_accessions(self, count: int) -> list[str]:
        """Return a list of pseudorandom, consecutive RefSeq accession numbers."""

        first_number = self.random_int(0, 999999 - count)

        accessions = [f"NC_{(first_number + i):06d}" for i in range(count)]

        return accessions


class BusinessProvider(BaseProvider):
    """Provider that generates business-related data like IDs."""

    def legacy_id(self) -> str:
        """Return a pseudorandom legacy ID."""
        return self.bothify("?????????")


class SequenceProvider(BaseProvider):
    """Dummy sequence data provider."""

    def sequence(
        self,
        min: int = MIN_SEQUENCE_LENGTH,
        max: int = MAX_SEQUENCE_LENGTH,
    ) -> str:
        """Return a pseudorandom string consisting of
        acceptable genetic sequence letters.
        """
        return "".join(
            self.random_elements(
                NUCLEOTIDE_PROBABILITIES,
                self.random_int(min, max),
                use_weighting=True,
            ),
        )

    def sequence_length(
        self,
        min: int = MIN_SEQUENCE_LENGTH,
        max: int = MAX_SEQUENCE_LENGTH,
    ) -> int:
        return self.random_int(min, max)


class OrganismProvider(BaseProvider):
    """Organism name provider. Recombines parts of preexisting taxon names
    to create quasi-realistic organisms.
    """

    def condition_adjective(self) -> str:
        """Return an adjective used to describe properties of an organism."""
        return self.random_element(ORGANISM_DESCRIPTOR_ADJECTIVES)

    def condition_noun(self) -> str:
        """Return a noun used to describe properties of an organism."""
        return self.random_element(ORGANISM_DESCRIPTOR_NOUNS)

    def host(self) -> str:
        """Return the host affected by the organism."""
        return self.random_element(ORGANISM_HOSTS)

    def part(self) -> str:
        """Return the part of the host affected by the organism."""
        return self.random_element(ORGANISM_PARTS)

    def virus_species(self) -> str:
        """Return the species of virus."""
        return self.random_element(ORGANISM_VIRUSES)

    def virus_type(self) -> str:
        """Return the type of organism."""
        return self.random_element(ORGANISM_TYPES)

    def host_part_and_descriptor_virus_organism(self) -> str:
        """Return an organism name consisting of HOST PART_AND_DESCRIPTOR VIRUS."""
        return " ".join(
            [
                self.host(),
                self.random_element(ORGANISM_PART_AND_DESCRIPTORS),
                self.virus_type(),
            ],
        ).capitalize()

    def host_adjective_virus_organism(self) -> str:
        """Return an organism name consisting of HOST ADJECTIVE VIRUS."""
        return (
            f"{self.host()} {self.condition_adjective()} {self.virus_type()}"
        ).capitalize()

    def host_noun_virus_organism(self):
        """Return an organism name consisting of HOST NOUN VIRUS."""
        return (
            f"{self.host()} {self.condition_noun()} {self.virus_type()}"
        ).capitalize()

    def host_adjective_noun_virus_organism(self):
        """Return an organism name consisting of HOST ADJECTIVE NOUN VIRUS."""
        return (
            f"{self.host()} {self.condition_adjective()} {self.condition_noun()} "
            f"{self.virus_type()}"
        ).capitalize()

    def organism(self):
        """Return a pseudorandom organism name."""
        match self.random_int(0, 9):
            case 0:
                return self.host_adjective_virus_organism()
            case 1:
                return self.host_adjective_noun_virus_organism()
            case 2:
                return self.host_noun_virus_organism()

        return self.host_part_and_descriptor_virus_organism()


class SegmentProvider(BaseProvider):
    """Quasi-realistic NCBISource segment field provider."""

    def segment_delimiter(self) -> str:
        """Return a segment delimiter."""
        return self.random_element({" ", "-", "_"})

    def segment_prefix(self) -> str:
        """Return a segment prefix."""
        return self.random_element({"DNA", "RNA"})

    def segment_key(self) -> str:
        """Return a segment key denoting the key identifier
        in the segment.
        """
        return self.random_element(
            {"1", "2", "A", "B", "C", "L", "N", "M", "R", "S", "U3"},
        )

    def segment(self) -> str:
        """Return a segment name."""
        if self.random_int(0, 9) < 7:
            return self.segment_prefix() + self.segment_delimiter() + self.segment_key()

        return self.segment_key()
