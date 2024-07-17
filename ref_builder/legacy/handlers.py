import json

from ref_builder.legacy.utils import (
    ErrorHandledResult,
    HandleErrorContext,
    extract_isolate_source,
)


def handle_enum(ctx: HandleErrorContext) -> ErrorHandledResult:
    """Handle a Pydantic ``enum`` error.

    If the field is ``source_type`` and fixing is enabled, the accessions associated
    with the isolate will be fetched from NCBI and the source name and source type will
    be derived and updated in the isolate.

    If no source type can be found or derived, the source_type will be set to "genbank"
    and the source_name will be set to the earliest sorted accession for the associated
    sequences.

    :param ctx: the error context
    :return: the error handling result

    """
    expected = ctx.error["ctx"]["expected"]
    field = ctx.error["loc"][-1]
    value = ctx.error["input"]

    message = (
        f"[bold]{field}[/bold] must be one of [u]{expected}[/u], but is "
        f"[u]{json.dumps(value).replace('"', "'")}[/u]."
    )

    if field == "source_type" and ctx.fix:
        isolate_index = ctx.error["loc"][1]
        isolate = ctx.otu["isolates"][isolate_index]

        genbank_records = ctx.ncbi_client.fetch_genbank_records(
            [s["accession"] for s in isolate["sequences"]],
        )

        try:
            source = extract_isolate_source(genbank_records)
        except ValueError:
            return ErrorHandledResult(
                f"{message}. no isolate source could be found on genbank.",
                False,
            )

        ctx.update_isolate(
            isolate_index,
            {"source_name": source.name, "source_type": source.type.lower()},
        )

        return ErrorHandledResult(message, True)

    return ErrorHandledResult(message)


def handle_int_type(ctx: HandleErrorContext) -> ErrorHandledResult:
    """Handle a Pydantic ``int_type`` error.

    If the field is ``taxid`` and fixing is enabled, the taxid will be fetched from NCBI
    using the OTU name and updated in the OTU.

    :param ctx: the error context
    :return: the error handling result
    """
    field = ctx.error["loc"][-1]
    value = ctx.error["input"]

    result = ErrorHandledResult(
        f"[bold]{field}[/bold] must be an integer. invalid value is "
        f"[u]{json.dumps(value)}[/u].",
    )

    if field == "taxid" and ctx.fix:
        taxid = ctx.ncbi_client.fetch_taxonomy_id_by_name(ctx.otu["name"])

        if taxid is None:
            return ErrorHandledResult(
                f"[bold]taxid[/bold] was not an integer and a taxid could not be "
                f"found for [u]{ctx.otu['name']}[/u] on ncbi.",
            )

        ctx.update_otu({"taxid": int(taxid)})

        result.fixed = True

        return result

    return result


def handle_missing(ctx: HandleErrorContext) -> ErrorHandledResult:
    """Handle a Pydantic ``missing`` error.

    If the field is ``taxid`` and fixing is enabled, the taxid will be fetched from NCBI
    using the OTU name and updated in the OTU.

    :param ctx: the error context
    :return: the error handling result

    """
    field = ctx.error["loc"][-1]

    message = f"[bold]{field}[/bold] is missing"

    if field == "taxid" and ctx.fix:
        taxid = ctx.ncbi_client.fetch_taxonomy_id_by_name(ctx.otu["name"])

        if taxid is None:
            return ErrorHandledResult(
                f"{message} and a taxid could not be found for [u]{ctx.otu['name']}[/u]"
                f" on ncbi.",
            )

        ctx.update_otu({"taxid": taxid})

        return ErrorHandledResult(message, fixed=True)

    return ErrorHandledResult(message)


def handle_str_type(ctx: HandleErrorContext) -> ErrorHandledResult:
    """Handle a Pydantic ``str_type`` error.

    If the field is ``accession`` and fixing is enabled, the accession will be fetched
    from NCBI using the OTU name and updated in the OTU.

    If the field is ``host`` and fixing is enabled, the host will be set to an empty
    string.

    :param ctx: the error context
    :return: the error handling result

    """
    value = ctx.error["input"]
    field = ctx.error["loc"][-1]

    if field == "host":
        return ErrorHandledResult(
            f"[bold]{field}[/bold] must be an string. invalid value is "
            f'[u]{json.dumps(value)}[/u]. set to [u]""[/u].',
            fixed=True,
        )

    return ErrorHandledResult(
        f"[bold]{field}[/bold] must be an string. invalid value is "
        f"[u]{json.dumps(value)}[/u].",
    )


def handle_string_pattern_mismatch(ctx: HandleErrorContext) -> ErrorHandledResult:
    """Handle a Pydantic ``string_pattern_mismatch`` error.

    These errors occur when an accession or genomic sequence do not match defined regex
    patterns.

    If the field is ``accession`` and the value does not contain a version suffix, the
    accession will be fetched from NCBI and the versioned accession will be attached
    to the sequence. Any other fields associated with the sequence will also be updated.

    If the field is ``sequence``, the error message will be updated to include the
    incorrect sequence characters.

    :param ctx: the error context
    :return: the error handling result
    """
    field = ctx.error["loc"][-1]
    pattern = ctx.error["ctx"]["pattern"]
    value = ctx.error["input"]

    if field == "sequence":
        invalid_chars = [
            (i, char) for i, char in enumerate(value) if char not in "ATCGRYKMSWBDHVN"
        ]

        mismatches = ",".join([f"{char}@{i}" for i, char in invalid_chars])

        return ErrorHandledResult(
            f"[bold]{field}[/bold] contains invalid characters. invalid "
            f"characters are [u]{mismatches}[/u].",
        )

    if field == "accession" and "." not in value:
        message = (
            f"[bold]accession[/bold] value [u]{value}[/u] does not contain a "
            f"[u].[/u] and needs a valid version suffix."
        )

        if not ctx.fix:
            return ErrorHandledResult(message)

        genbank_records = ctx.ncbi_client.fetch_genbank_records(
            [value],
        )

        if not genbank_records:
            return ErrorHandledResult(
                message + " couldn't fix because no matching genbank record found.",
            )

        genbank_record = genbank_records[0]

        isolate_index = ctx.error["loc"][1]
        sequence_index = ctx.error["loc"][3]

        ctx.update_sequence(
            isolate_index,
            sequence_index,
            {
                "accession": genbank_record.accession_version,
                "definition": genbank_record.definition,
                "sequence": genbank_record.sequence,
            },
        )

        return ErrorHandledResult(message, fixed=True)

    return ErrorHandledResult(
        f"[bold]{field}[/bold] does not match the required regex pattern "
        f"[u]{pattern}[/u].",
    )


def handle_string_too_short(ctx: HandleErrorContext) -> ErrorHandledResult:
    """Handle a Pydantic ``string_too_short`` error.

    If the field is ``source_name`` and fixing is enabled, the accessions associated
    with the isolate will be fetched from NCBI and the source name and source type will
    be derived and updated in the isolate.

    If no source type can be found or derived, an informative error message will be
    returned.

    :param ctx: the error context
    :return: the error handling result
    """
    field = ctx.error["loc"][-1]
    input_ = ctx.error["input"]
    min_length = ctx.error["ctx"]["min_length"]

    result = ErrorHandledResult(
        f"[bold]{field}[/bold] must be a string with a minimum length of "
        f"{min_length}. it currently has a length of [u]{len(input_)}[/u].",
    )

    if field == "source_name" and ctx.fix:
        isolate_index = ctx.error["loc"][1]
        isolate = ctx.otu["isolates"][isolate_index]

        genbank_records = ctx.ncbi_client.fetch_genbank_records(
            [s["accession"] for s in isolate["sequences"]],
        )

        if not genbank_records:
            return ErrorHandledResult(
                f"[bold]source_name[/bold] must be a string with a minimum length of "
                f"{min_length}. it currently has a length of [u]{len(input_)}[/u]. we "
                "could not find a source name on genbank.",
            )

        source = extract_isolate_source(genbank_records)

        ctx.update_isolate(
            isolate_index,
            {"source_name": source.name, "source_type": source.type.lower()},
        )

        result.fixed = True

    return result


def handle_too_short(ctx: HandleErrorContext) -> ErrorHandledResult:
    """Handle a Pydantic ``too_short`` error.

    Return a result object with an informative, human-readable error message.

    :param ctx: the error context
    :return: the error handling result
    """
    actual_length = ctx.error["ctx"]["actual_length"]
    field = ctx.error["loc"][-1]
    field_type = ctx.error["ctx"]["field_type"]
    min_length = ctx.error["ctx"]["min_length"]

    return ErrorHandledResult(
        f"[bold]{field}[/bold] must be a {field_type} with a minimum length of "
        f"{min_length}. it currently has a length of [u]{actual_length}[/u].",
    )


def handle_value_error(ctx: HandleErrorContext) -> ErrorHandledResult:
    """Handle a ``ValueError`` by returning a result object with a human-readable error
    message.

    :param ctx: the error context
    :return: the error handling result
    """
    return ErrorHandledResult(str(ctx.error["ctx"]["error"]))
