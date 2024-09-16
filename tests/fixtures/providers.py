from collections import OrderedDict

from faker.providers import BaseProvider, DynamicProvider

from ref_builder.utils import Accession
from tests.fixtures.data.organism import (
    ORGANISM_HOSTS,
    ORGANISM_PARTS,
    ORGANISM_DESCRIPTOR_NOUNS,
    ORGANISM_DESCRIPTOR_ADJECTIVES,
    ORGANISM_PART_AND_DESCRIPTORS,
    ORGANISM_TYPES,
    ORGANISM_VIRUSES,
    ORGANISM_SPECIFIC_VIRUS_TYPES,
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
    def refseq_accession(self):
        return "NC_" + self.numerify("######")

    def genbank_accession(self):
        dice_roll = self.random_int(0, 10)

        if dice_roll > 6:
            return self.bothify("?#####").upper()

        return self.bothify("??######").upper()

    def accession(self) -> str:
        dice_roll = self.random_int(0, 10)

        if dice_roll > 7:
            return self.genbank_accession()
        else:
            return self.refseq_accession()


class SequenceProvider(BaseProvider):
    def sequence(self):
        return "".join(
            self.random_elements(
                SEQUENCE_DICTIONARY,
                self.random_int(SEQUENCE_MIN, SEQUENCE_MAX),
                use_weighting=True,
        ))


class OrganismProvider(BaseProvider):
    def __init__(self, generator):
        super().__init__(generator)

    def condition_adjective(self):
        return self.random_element(ORGANISM_DESCRIPTOR_ADJECTIVES)

    def condition_noun(self):
        return self.random_element(ORGANISM_DESCRIPTOR_NOUNS)

    def host(self):
        return self.random_element(ORGANISM_HOSTS)

    def part(self):
        return self.random_element(ORGANISM_PARTS)

    def virus_species(self):
        return self.random_element(ORGANISM_VIRUSES)

    def virus_type(self):
        return self.random_element(ORGANISM_TYPES)

    def host_part_and_descriptor_virus_organism(self):
        return " ".join(
            [
                self.host(),
                self.random_element(ORGANISM_PART_AND_DESCRIPTORS),
                self.virus_type(),
            ]
        ).capitalize()

    def host_adjective_virus_organism(self):
        return " ".join(
            [
                self.host(),
                self.condition_adjective(),
                self.virus_type(),
            ]
        ).capitalize()

    def host_noun_virus_organism(self):
        return " ".join(
            [
                self.host(),
                self.condition_noun(),
                self.virus_type(),
            ]
        ).capitalize()

    def host_adjective_noun_virus_organism(self):
        return " ".join(
            [
                self.host(),
                self.condition_adjective(),
                self.condition_noun(),
                self.virus_type(),
            ]
        ).capitalize()

    def organism(self):
        roll = self.random_int(0, 9)

        match roll:
            case 0:
                return self.host_adjective_virus_organism()
            case 1:
                return self.host_adjective_noun_virus_organism()
            case 2:
                return self.host_noun_virus_organism()

        return self.host_part_and_descriptor_virus_organism()


class SourceProvider(BaseProvider):
    def __init__(self, generator):
        super().__init__(generator)

    def segment_delimiter(self):
        return self.random_element({" ", "-", "_"})

    def segment_prefix(self):
        return self.random_element({"DNA", "RNA"})

    def segment_key(self):
        return self.random_element({"1", "2", "A", "B", "C", "L", "N", "M", "R", "S", "U3"})

    def segment(self):
        if self.random_int(0, 9) < 7:
            return self.segment_prefix() + self.segment_delimiter() + self.segment_key()

        return self.segment_key()


if __name__ == '__main__':
    from faker import Faker

    fake = Faker()
    fake.add_provider(SourceProvider)
    fake.add_provider(OrganismProvider)
    fake.add_provider(SequenceProvider)
    fake.add_provider(AccessionProvider)

    for _ in range(10):
        print(fake.segment())
        print(fake.organism())
        print(fake.sequence())
        print(fake.accession())
        print()

