"""A client for fetching data from NCBI databases."""

import os
from collections.abc import Collection
from contextlib import contextmanager
import datetime
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

base_logger = get_logger("ncbi")

ESEARCH_PAGE_SIZE = 1000
"""The number of results to fetch per page in an Entrez esearch query."""

DATE_TEMPLATE = "%Y/%m/%d"
"""The standard date format used by NCBI Entrez."""


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

    def fetch_genbank_records(self, accessions: Collection[str]) -> list[NCBIGenbank]:
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
            uncached_accessions = []

            for accession in accessions:
                record = self.cache.load_genbank_record(accession)
                if record is not None:
                    records.append(record)

            if records:
                logger.debug(
                    f"Loaded {len(records)} cached records",
                    cached_accessions=[
                        record.get(GenbankRecordKey.PRIMARY_ACCESSION)
                        for record in records
                    ],
                )
            if uncached_accessions:
                logger.debug(
                    "Uncached accessions found",
                    uncached_accessions=uncached_accessions,
                )

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
    def fetch_unvalidated_genbank_records(accessions: Collection[str]) -> list[dict]:
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
                        id=list(accessions),
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
    def fetch_accessions_by_taxid(
        taxid: int,
        sequence_min_length: int = 0,
        sequence_max_length: int = 0,
        modification_date_start: datetime.date | None = None,
        modification_date_end: datetime.date | None = None,
        refseq_only: bool = False,
    ) -> list[str]:
        """Fetch all accessions associated with the given ``taxid``.

        :param taxid: A Taxonomy ID
        :param sequence_min_length: The minimum length of a fetched sequence.
        :param sequence_max_length: The maximum length of a fetched sequence.
        :param modification_date_start: The earliest a sequence's latest modification date can be.::
        :param modification_date_end: The latest a sequence's latest modification date can be.::
        :return: A list of Genbank accessions
        """
        page = 1
        accessions = []

        term = f"txid{taxid}[orgn]"
        if sequence_min_length > 0 and sequence_max_length > 0:
            term += " AND " + NCBIClient.generate_sequence_length_filter_string(
                sequence_min_length,
                sequence_max_length,
            )

        if modification_date_start is not None:
            term += " AND " + NCBIClient.generate_date_filter_string(
                "MDAT",
                modification_date_start,
                modification_date_end,
            )

        if refseq_only:
            term += " AND " + "refseq[filter]"

        # If there are more than 1000 accessions, we need to paginate.
        while True:
            retstart = (page - 1) * ESEARCH_PAGE_SIZE

            with log_http_error():
                handle = Entrez.esearch(
                    db=NCBIDatabase.NUCCORE,
                    term=term,
                    idtype="acc",
                    retstart=retstart,
                    retmax=ESEARCH_PAGE_SIZE,
                )

            result = Entrez.read(handle)

            if not result["IdList"]:
                break

            result_count = int(result["Count"])

            accessions += result["IdList"]

            if result_count - retstart <= ESEARCH_PAGE_SIZE:
                break

            base_logger.debug(
                "Large fetch. May take longer than expected...",
                result_count=result_count,
                page=page,
                page_size=ESEARCH_PAGE_SIZE,
                taxid=taxid,
            )

            page += 1

        return accessions

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

        try:
            return NCBITaxonomy.model_validate(record)

        except ValidationError:
            rank = self._fetch_taxonomy_rank(taxid)
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

    @staticmethod
    def filter_accessions(raw_accessions: Collection[str]) -> set[Accession]:
        """Filter raw eSearch accession list and return a set of compliant Nucleotide accessions."""
        valid_accessions = set()
        for accession in raw_accessions:
            try:
                valid_accessions.add(Accession.from_string(accession))
            except ValueError:
                pass

        return valid_accessions

    @staticmethod
    def generate_sequence_length_filter_string(
        min_length: int = 0, max_length: int = 0
    ) -> str:
        """Return a term filter string delimiting a given length range.

        Returns an empty string if not given a min or max length parameter.
        """
        if min_length > 0:
            if min_length > 0:
                return f'"{min_length}"[SLEN] : "{max_length}"[SLEN]'

            return f'"{min_length}"[SLEN]'

        if max_length > 0:
            return f'"{max_length}"[SLEN]'

        return ""

    @staticmethod
    def generate_date_filter_string(
        filter_type: str,
        start_date: datetime.date | None = None,
        end_date: datetime.date | None = None,
    ) -> str:
        """Return a term filter string delimiting a given time range.

        Returns an empty string if not given a start or end date parameter.

        :param filter_type: The search term filter type. Can be "MDAT" or "PDAT".:
        :param start_date: The start date of the search range. If None, assume no lower bound.:
        :param end_date: The end date of the search range. If None, assume no upper bound.:
        :return: A formatted date filter string or an empty string.
        """
        if filter_type not in {"MDAT", "PDAT"}:
            raise ValueError(
                "Invalid filter type. Only ``MDAT``, ``PDAT`` are supported."
            )

        if start_date is None and end_date is None:
            return ""

        start_date_string = "0001/01/01"
        end_date_string = "3000/12/31"

        if start_date:
            start_date_string = start_date.strftime(DATE_TEMPLATE)

            if end_date:
                end_date_string = end_date.strftime(DATE_TEMPLATE)

        return (
            f'"{start_date_string}"[{filter_type}]'
            + " : "
            + f'"{end_date_string}"[{filter_type}]'
        )


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
