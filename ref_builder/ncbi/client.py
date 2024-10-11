"""A client for fetching data from NCBI databases."""

import os
from contextlib import contextmanager
from enum import StrEnum
from http import HTTPStatus
from urllib.error import HTTPError
from urllib.parse import quote_plus

from Bio import Entrez
from pydantic import ValidationError
from structlog import get_logger

from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.models import NCBIDatabase, NCBIGenbank, NCBIRank, NCBITaxonomy
from ref_builder.utils import Accession

if email := os.environ.get("NCBI_EMAIL"):
    Entrez.email = email

if api_key := os.environ.get("NCBI_API_KEY"):
    Entrez.api_key = os.environ.get("NCBI_API_KEY")

base_logger = get_logger()

ESEARCH_PAGE_SIZE = 1000
"""The number of results to fetch per page in an Entrez esearch query."""


class GenbankRecordKey(StrEnum):
    """Genbank record keys."""

    ACCESSION_VERSION = "GBSeq_accession-version"
    COMMENT = "GBSeq_comment"
    DEFINITION = "GBSeq_definition"
    FEATURE_TABLE = "GBSeq_feature-table"
    LENGTH = "GBSeq_length"
    PRIMARY_ACCESSION = "GBSeq_primary-accession"
    SEQUENCE = "GBSeq_sequence"


class NCBIClient:
    """A client for fetching, caching, and validating data from NCBI databases."""

    def __init__(self, ignore_cache: bool) -> None:
        """Initialize the NCBI client with a cache path and an ignore_cache flag.

        :param ignore_cache: If True, does not allow the return of cached data
        """
        self.cache = NCBICache()
        self.ignore_cache = ignore_cache

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
            set(accessions)
            - {record.get(GenbankRecordKey.PRIMARY_ACCESSION) for record in records},
        )

        if fetch_list:
            logger.debug("Fetching accessions...", fetch_list=fetch_list)
            new_records = NCBIClient.fetch_unvalidated_genbank_records(fetch_list)

            for record in new_records:
                versioned_accession = Accession.from_string(
                    record[GenbankRecordKey.ACCESSION_VERSION],
                )
                self.cache.cache_genbank_record(
                    record,
                    versioned_accession.key,
                    versioned_accession.version,
                )

            if new_records:
                records.extend(new_records)

        if records:
            return sorted(
                NCBIClient.validate_genbank_records(records),
                key=lambda r: r.accession,
            )

        return []

    @staticmethod
    def fetch_unvalidated_genbank_records(accessions: list[str]) -> list[dict]:
        """Fetch a list of Genbank records given a list of accessions.

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
                    logger.warning("Bad ID.", exception=e)
                    return []

        except HTTPError as e:
            if e.code == HTTPStatus.BAD_REQUEST:
                logger.exception("Accessions not found")
            else:
                logger.exception()

            return []

        try:
            records = Entrez.read(handle)
        except RuntimeError:
            logger.exception("NCBI returned unparseable data")
            return []

        if records:
            # Handle cases where not all accessions can be fetched
            if len(records) != len(accessions):
                logger.debug("Partial results fetched. Returning results...")

            return records

        return []

    @staticmethod
    def fetch_accessions_by_taxid(taxid: int) -> list[str]:
        """Fetch all accessions associated with the given ``taxid``.

        :param taxid: A Taxonomy ID
        :return: A list of Genbank accessions
        """
        page = 1

        accessions = []

        # If there are more than 1000 accessions, we need to paginate.
        while True:
            with log_http_error():
                handle = Entrez.esearch(
                    db=NCBIDatabase.NUCCORE,
                    term=f"txid{taxid}[orgn]",
                    idtype="acc",
                    retstart=page * ESEARCH_PAGE_SIZE,
                    retmax=ESEARCH_PAGE_SIZE,
                )

            result = Entrez.read(handle)

            accessions.extend(result["IdList"])

            if int(result["Count"]) < ESEARCH_PAGE_SIZE:
                break

            page += 1

        with log_http_error():
            handle = Entrez.esearch(
                db=NCBIDatabase.NUCCORE,
                term=f"txid{taxid}[orgn]",
                idtype="acc",
                retstart=0,
                retmax=1000,
            )

        return Entrez.read(handle)["IdList"]

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
            accession = record.get(GenbankRecordKey.PRIMARY_ACCESSION, "?")

            try:
                clean_records.append(NCBIGenbank.model_validate(record))

            except ValidationError as exc:
                base_logger.debug(
                    "Encountered validation errors",
                    accession=accession,
                    count=exc.error_count(),
                    errors=exc.errors(),
                )
                continue

        return clean_records

    def fetch_taxonomy_record(self, taxid: int) -> NCBITaxonomy | None:
        """Fetch and validate a taxonomy record from NCBI Taxonomy.

        If the record rank has an invalid rank (e.g. "no data"), makes an additional
        docsum fetch and attempts to extract the rank data.

        :param taxid: A NCBI Taxonomy id
        :return: A validated NCBI Taxonomy record NCBITaxonomy if possible,
            else None
        """
        logger = base_logger.bind(taxid=taxid)

        record = None if self.ignore_cache else self.cache.load_taxonomy(taxid)

        if record is None:
            with log_http_error():
                try:
                    handle = Entrez.efetch(
                        db=NCBIDatabase.TAXONOMY,
                        id=taxid,
                        rettype="null",
                    )
                except HTTPError:
                    return None

            try:
                records = Entrez.read(handle)
            except RuntimeError:
                logger.exception("NCBI returned unparseable data.")
                return None

            if records:
                record = records[0]
                self.cache.cache_taxonomy_record(record, taxid)
            else:
                return None

        if rank := self._fetch_taxonomy_rank(taxid):
            try:
                return NCBITaxonomy(**record, rank=rank)
            except ValidationError:
                logger.exception("Failed to find a valid rank.")

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
        except HTTPError:
            return None

        try:
            docsum_record = Entrez.read(handle)
        except RuntimeError:
            logger.exception("NCBI returned unparseable data")
            return None

        try:
            rank = NCBIRank(docsum_record[0]["Rank"])
        except ValueError:
            logger.debug("Found rank for this taxid, but it did not pass validation.")
            return None

        logger.debug("Found valid rank", rank=rank)

        return rank

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
        except HTTPError:
            return None

        try:
            record = Entrez.read(handle)
        except RuntimeError:
            logger.exception("NCBI returned unparseable data.")
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
        except HTTPError:
            return None

        try:
            record = Entrez.read(handle)
        except RuntimeError:
            logger.exception("NCBI returned unparseable data.")
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
        base_logger.exception(
            "HTTPError raised",
            body=e.read(),
            code=e.code,
            reason=e.reason,
        )
        raise
