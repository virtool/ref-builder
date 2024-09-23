from collections import OrderedDict

from faker.providers import BaseProvider

from tests.fixtures.data.organism import (
    ORGANISM_DESCRIPTOR_ADJECTIVES,
    ORGANISM_DESCRIPTOR_NOUNS,
    ORGANISM_HOSTS,
    ORGANISM_PART_AND_DESCRIPTORS,
    ORGANISM_PARTS,
    ORGANISM_TYPES,
    ORGANISM_VIRUSES,
)

UNCOMMON_PROBABILITY = 0.01
SEQUENCE_DICTIONARY = OrderedDict([
    ("A", 0.245),
    ("T", 0.245),
    ("C", 0.245),
    ("G", 0.245),
    ("R", UNCOMMON_PROBABILITY),
    ("Y", UNCOMMON_PROBABILITY),
    ("K", UNCOMMON_PROBABILITY),
    ("M", UNCOMMON_PROBABILITY),
    ("S", UNCOMMON_PROBABILITY),
    ("W", UNCOMMON_PROBABILITY),
    ("B", UNCOMMON_PROBABILITY),
    ("D", UNCOMMON_PROBABILITY),
    ("H", UNCOMMON_PROBABILITY),
    ("V", UNCOMMON_PROBABILITY),
    ("N", UNCOMMON_PROBABILITY),
])
SEQUENCE_MIN = 100
SEQUENCE_MAX = 1500


class AccessionProvider(BaseProvider):
    """Raw accession provider based on GenBank's guidelines for accession numbers."""

    def refseq_accession(self) -> str:
        """Return a pseudorandom RefSeq accession number."""
        return "NC_" + self.numerify("######")

    def genbank_accession(self) -> str:
        """Return a pseudorandom non-RefSeq accession number."""
        dice_roll = self.random_int(0, 10)

        if dice_roll > 6:
            return self.bothify("?#####").upper()

        return self.bothify("??######").upper()

    def accession(self) -> str:
        """Return a pseudorandom accession number."""
        dice_roll = self.random_int(0, 10)

        if dice_roll > 7:
            return self.genbank_accession()

        return self.refseq_accession()


class SequenceProvider(BaseProvider):
    """Dummy sequence data provider."""

    def sequence(self) -> str:
        """Return a pseudorandom string consisting of
        acceptable genetic sequence letters.
        """
        return "".join(
            self.random_elements(
                SEQUENCE_DICTIONARY,
                self.random_int(SEQUENCE_MIN, SEQUENCE_MAX),
                use_weighting=True,
        ))


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
            ]
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
            f"{self.host()} {self.condition_adjective()} {self.condition_noun()} {self.virus_type()}"
        ).capitalize()

    def organism(self):
        """Return a pseudorandom organism name."""
        roll = self.random_int(0, 9)

        match roll:
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
