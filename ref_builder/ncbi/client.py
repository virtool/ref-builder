import os
from contextlib import contextmanager
from enum import StrEnum
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import quote_plus

from Bio import Entrez
from pydantic import ValidationError
from structlog import get_logger

from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.models import NCBIDatabase, NCBIGenbank, NCBIRank, NCBITaxonomy

if email := os.environ.get("NCBI_EMAIL"):
    Entrez.email = email
if api_key := os.environ.get("NCBI_API_KEY"):
    Entrez.api_key = os.environ.get("NCBI_API_KEY")

base_logger = get_logger()


class GBSeq(StrEnum):
    ACCESSION = "GBSeq_primary-accession"
    DEFINITION = "GBSeq_definition"
    SEQUENCE = "GBSeq_sequence"
    LENGTH = "GBSeq_length"
    COMMENT = "GBSeq_comment"
    FEATURE_TABLE = "GBSeq_feature-table"


class NCBIClient:
    def __init__(self, cache_path: Path, ignore_cache: bool) -> None:
        """Initialize the NCBI client with a cache path and an ignore_cache flag.

        :param cache_path: A path to a directory to be used as a cache
        :param ignore_cache: If True, does not allow the return of cached data
        """
        self.cache = NCBICache(cache_path)
        self.ignore_cache = ignore_cache

    @classmethod
    def from_repo(cls, path: Path, ignore_cache: bool) -> "NCBIClient":
        """Initialize the NCBI cache in the default subpath under a given repository.

        :param path: the path to a reference repository
        :param ignore_cache: whether cached data should be ignored
        :return: an instance of NCBIClient
        """
        return NCBIClient(path / ".cache/ncbi", ignore_cache=ignore_cache)

    def fetch_genbank_records(self, accessions: list[str]) -> list[NCBIGenbank]:
        """Fetch or load NCBI Nucleotide records corresponding to a list of accessions.

        Cache fetched records if found. Returns validated records.

        :param accessions: A list of accessions to be fetched
        :return: A list of validated NCBIGenbank records
        """
        if not accessions:
            return []

        records = []

        logger = base_logger.bind(accessions=accessions)

        if not self.ignore_cache:
            for accession in accessions:
                record = self.cache.load_genbank_record(accession)
                if record is not None:
                    records.append(record)
                else:
                    logger.debug("Missing accession", missing_accession=accession)

        fetch_list = list(
            set(accessions) - {record.get(GBSeq.ACCESSION) for record in records},
        )

        if fetch_list:
            logger.debug("Fetching accessions...", fetch_list=fetch_list)
            new_records = NCBIClient.fetch_unvalidated_genbank_records(fetch_list)

            for record in new_records:
                try:
                    self.cache.cache_genbank_record(record, record[GBSeq.ACCESSION])
                except FileNotFoundError:
                    logger.error("Failed to cache record")

            if new_records:
                records.extend(new_records)

        if records:
            return NCBIClient.validate_genbank_records(records)

        return []

    @staticmethod
    def fetch_unvalidated_genbank_records(accessions: list[str]) -> list[dict]:
        """Take a list of accession numbers, parse the corresponding XML records
        from GenBank using Entrez.Parser and return

        :param accessions: List of accession numbers to fetch from GenBank
        :return: A list of deserialized XML records from NCBI Nucleotide
        """
        logger = base_logger.bind(accessions=accessions)

        try:
            with log_http_error():
                try:
                    handle = Entrez.efetch(
                        db=NCBIDatabase.NUCCORE,
                        id=accessions,
                        rettype="gb",
                        retmode="xml",
                    )
                except RuntimeError as e:
                    logger.warning(f"Bad ID. (RunTime Error: {e})")
                    return []

        except HTTPError as e:
            if e.code == 400:
                logger.error(f"Accessions not found ({e})")
            else:
                logger.error(e)

            return []

        try:
            records = Entrez.read(handle)
        except RuntimeError as e:
            logger.error(f"NCBI returned unparseable data ({e})")
            return []

        if records:
            # Handle cases where not all accessions can be fetched
            if len(records) != len(accessions):
                logger.debug("Partial results fetched. Returning results...")

            return records

        return []

    def link_from_taxid_and_fetch(self, taxid: int) -> list[NCBIGenbank]:
        """Fetch all Genbank records linked to a taxonomy record.

        Usable without preexisting OTU data.

        :param taxid: A NCBI Taxonomy id
        :return: A list of Entrez-parsed Genbank records
        """
        accessions = NCBIClient.link_accessions_from_taxid(taxid)

        logger = base_logger.bind(taxid=taxid, linked_accessions=accessions)

        logger.debug("Fetching accessions...")
        records = NCBIClient.fetch_unvalidated_genbank_records(accessions)

        if not records:
            logger.info("No accessions found for this Taxon ID")
            return []

        for record in records:
            try:
                self.cache.cache_genbank_record(record, record[GBSeq.ACCESSION])
            except FileNotFoundError:
                logger.error("Failed to cache record")

        if records:
            return NCBIClient.validate_genbank_records(records)

        return []

    @staticmethod
    def link_accessions_from_taxid(taxid: int) -> list:
        """Requests a cross-reference for NCBI Taxonomy and Nucleotide via ELink
        and returns the results as a list.

        :param taxid: A NCBI Taxonomy id
        :return: A list of Genbank accessions linked to the Taxonomy UID
        """
        logger = base_logger.bind(taxid=taxid)
        try:
            with log_http_error():
                handle = Entrez.elink(
                    dbfrom="taxonomy",
                    db=NCBIDatabase.NUCCORE,
                    id=str(taxid),
                    idtype="acc",
                )
        except HTTPError as e:
            logger.error(f"{e.code}: {e.reason}")
            logger.error("Your request was likely refused by NCBI.")
            return []

        try:
            elink_results = Entrez.read(handle)
        except RuntimeError as e:
            logger.error(f"NCBI returned unparseable data ({e})")
            return []

        if elink_results:
            # Discards unneeded tables and formats needed table as a list
            for link_set_db in elink_results[0]["LinkSetDb"]:
                if link_set_db["LinkName"] == "taxonomy_nuccore":
                    id_table = link_set_db["Link"]

                    return [keypair["Id"].split(".")[0] for keypair in id_table]

        return []

    @staticmethod
    def validate_genbank_records(records: list[dict]) -> list[NCBIGenbank]:
        """Process a list of raw Genbank dicts into validated NCBIGenbank records.
        Logs an error if there is an issue with validation or parsing,
        but does not fail out.

        :param records: A list of unvalidated NCBI Genbank records
        :return: A list of validated records as NCBIGenbank
        """
        clean_records = []

        for record in records:
            accession = record.get(GBSeq.ACCESSION, "?")

            try:
                clean_records.append(NCBIClient.validate_genbank_record(record))

            except ValidationError as exc:
                base_logger.debug(
                    f"Encountered {exc.error_count()} validation errors",
                    accession=accession,
                    errors=exc.errors(),
                )
                continue

        return clean_records

    @staticmethod
    def validate_genbank_record(raw: dict) -> NCBIGenbank:
        """Parses an NCBI Genbank record from a Genbank dict to a validated NCBIGenbank.

        :param raw: A NCBI Nucleotide dict record, parsed by Bio.Entrez.Parser
        :return: A validated subset of Genbank record data
        """
        return NCBIGenbank(**raw)

    def fetch_taxonomy_record(self, taxid: int) -> NCBITaxonomy | None:
        """Fetch and validate a taxonomy record from NCBI Taxonomy.

        If the record rank has an invalid rank (e.g. "no data"), makes an additional
        docsum fetch and attempts to extract the rank data.

        :param taxid: A NCBI Taxonomy id
        :return: A validated NCBI Taxonomy record NCBITaxonomy if possible,
            else None
        """
        logger = base_logger.bind(taxid=taxid)

        record = None
        if not self.ignore_cache:
            record = self.cache.load_taxonomy(taxid)
            if record:
                logger.debug("Cached record found")
            else:
                logger.debug("Cached record not found")

        if record is None:
            with log_http_error():
                try:
                    handle = Entrez.efetch(
                        db=NCBIDatabase.TAXONOMY,
                        id=taxid,
                        rettype="null",
                    )
                except HTTPError as e:
                    logger.error(f"{e.code}: {e.reason}")
                    logger.error("Your request was likely refused by NCBI.")
                    return None

            try:
                records = Entrez.read(handle)
            except RuntimeError as e:
                logger.error(f"NCBI returned unparseable data ({e})")
                return None

            if records:
                record = records[0]
                logger.debug("Caching data...")
                self.cache.cache_taxonomy_record(record, taxid)
            else:
                logger.error("ID not found in NCBI Taxonomy database.")
                return None

        if record is None:
            return None

        try:
            return NCBIClient.validate_taxonomy_record(record)
        except ValidationError as exc:
            for error in exc.errors():
                if error["loc"][0] == "Rank":
                    logger.warning(
                        "Rank data not found in record",
                        input=error["input"],
                        loc=error["loc"][0],
                        msg=error["msg"],
                    )
                else:
                    logger.error(
                        "Taxonomy record failed validation",
                        errors=exc.errors(),
                    )
                    return None

        logger.info("Running additional docsum fetch...")

        try:
            rank = self._fetch_taxonomy_rank(taxid)
        except HTTPError as e:
            logger.error(f"{e.code}: {e.reason}")
            logger.error("Your request was likely refused by NCBI.")
            return None

        if rank:
            try:
                return NCBIClient.validate_taxonomy_record(record, rank)
            except ValidationError as exc:
                logger.error("Failed to find a valid rank.", error=exc)

        return None

    @staticmethod
    def _fetch_taxonomy_rank(taxid: int) -> NCBIRank | None:
        """Fetch the rank for a given NCBI Taxonomy ID.

        If no rank is found in the record, ``None`` is returned.

        :param taxid: an NCBI Taxonomy ID
        :return: the rank if possible
        """
        logger = base_logger.bind(taxid=taxid)

        try:
            with log_http_error():
                handle = Entrez.efetch(
                    db=NCBIDatabase.TAXONOMY,
                    id=taxid,
                    rettype="docsum",
                    retmode="xml",
                )
        except HTTPError as e:
            logger.error(f"{e.code}: {e.reason}")
            logger.error("Your request was likely refused by NCBI.")
            return None

        try:
            docsum_record = Entrez.read(handle)
        except RuntimeError as e:
            logger.error(f"NCBI returned unparseable data ({e})")
            return None

        try:
            rank = NCBIRank(docsum_record[0]["Rank"])
        except ValueError:
            logger.debug("Found rank for this taxid, but it did not pass validation")
            return None

        logger.debug("Valid rank found", rank=rank)
        return rank

    @staticmethod
    def validate_taxonomy_record(
        record: dict,
        rank: NCBIRank | None = None,
    ) -> NCBITaxonomy:
        """Attempt to validate a raw record from NCBI Taxonomy.

        If the ``rank`` parameter is set, it will override the inline rank data.

        Throws a ValidationError if valid rank data is not found in the record.

        :param record: Unvalidated taxonomy record data
        :param rank: Overrides inline rank data if set.
            Use to insert a valid rank from a docsum fetch
            on a second validation attempt.
        :return: A validated NCBITaxonomy record
        """
        if rank is None:
            return NCBITaxonomy(**record)

        return NCBITaxonomy(rank=rank, **record)

    @staticmethod
    def fetch_taxonomy_id_by_name(name: str) -> int | None:
        """Fetch a best-guess taxonomy ID for a given OTU name.

        Searches the NCBI taxonomy database for a given OTU name, then fetches and
        returns its taxonomy ID.

        Returns ``None`` if no matching taxonomy is found.

        :param name: The name of an otu
        :return: The NCBI Taxonomy UID for the given otu name
        """
        logger = base_logger.bind(name=name)
        try:
            with log_http_error():
                handle = Entrez.esearch(db="taxonomy", term=name)
        except HTTPError as e:
            logger.error(f"{e.code}: {e.reason}")
            logger.error("Your request was likely refused by NCBI.")
            return None

        try:
            record = Entrez.read(handle)
        except RuntimeError as e:
            logger.error(f"NCBI returned unparseable data ({e})")
            return None

        try:
            return int(record["IdList"][0])
        except IndexError:
            return None

    @staticmethod
    def fetch_spelling(
        name: str,
        db: NCBIDatabase = NCBIDatabase.TAXONOMY,
    ) -> str | None:
        """Fetch alternative spellings for a given OTU name.

        :param name: The OTU name that requires correcting
        :param db: Database to check against. Defaults to ``taxonomy``.
        :return: String containing NCBI-suggested spelling changes, None if HTTPError
        """
        logger = base_logger.bind(query=name)
        try:
            with log_http_error():
                handle = Entrez.espell(db=db, term=quote_plus(name))
        except HTTPError as e:
            logger.error(f"{e.code}: {e.reason}")
            logger.error("Your request was likely refused by NCBI.")
            return None

        try:
            record = Entrez.read(handle)
        except RuntimeError as e:
            logger.error(f"NCBI returned unparseable data ({e})")
            return None

        if "CorrectedQuery" in record:
            return record["CorrectedQuery"]

        logger.warning("Suggested spelling not found.")
        return None


@contextmanager
def log_http_error() -> None:
    """Log detailed HTTPError info for debugging before throwing the HTTPError."""
    try:
        yield

    except HTTPError as e:
        base_logger.error(
            "HTTPError raised",
            code=e.code,
            reason=e.reason,
            body=e.read(),
        )
        base_logger.exception(e)
        raise e
