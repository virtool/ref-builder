"""Validates a legacy repo.

A legacy repo must be validated before being converted to an event-sourced repo.
"""

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from pydantic import ValidationError
from rich.padding import Padding
from rich.table import Table
from structlog import get_logger

from ref_builder.console import console
from ref_builder.legacy.handlers import (
    handle_enum,
    handle_int_type,
    handle_missing,
    handle_str_type,
    handle_string_pattern_mismatch,
    handle_string_too_short,
    handle_too_short,
    handle_value_error,
)
from ref_builder.legacy.models import LegacyOTU
from ref_builder.legacy.repo import (
    check_unique_accessions,
    check_unique_ids,
    check_unique_otu_abbreviations_and_names,
)
from ref_builder.legacy.utils import (
    ErrorHandledResult,
    HandleErrorContext,
    build_legacy_otu,
    replace_otu,
)
from ref_builder.ncbi.client import NCBIClient

logger = get_logger("legacy")


@dataclass
class OTUValidationResult:
    handler_results: list[ErrorHandledResult]
    repaired_otu: dict


def handle_validation_error(
    e: ValidationError,
    fix: bool,
    ncbi_client: NCBIClient,
    otu: dict,
) -> OTUValidationResult:
    """Turn a ValidationError into a human-readable error message.

    :param e: the validation error
    :param fix: whether to attempt to fix the errors
    :param ncbi_client: the NCBI client
    :param otu: the otu data
    """
    handler_results = []
    repaired_otu = deepcopy(otu)

    for error in e.errors():
        error_type = error["type"]

        try:
            handler = {
                "enum": handle_enum,
                "int_type": handle_int_type,
                "missing": handle_missing,
                "string_type": handle_str_type,
                "string_pattern_mismatch": handle_string_pattern_mismatch,
                "string_too_short": handle_string_too_short,
                "too_short": handle_too_short,
                "value_error": handle_value_error,
            }[error_type]
        except KeyError:
            result = ErrorHandledResult(
                f"unhandled error type: {error_type}",
                False,
            )
        else:
            result = handler(
                HandleErrorContext(error, fix, ncbi_client, otu, repaired_otu),
            )

        handler_results.append(result)

    return OTUValidationResult(
        handler_results=handler_results,
        repaired_otu=repaired_otu,
    )


def validate_legacy_otu(
    fix: bool,
    ncbi_client: NCBIClient,
    otu: dict,
) -> OTUValidationResult | None:
    """Validate a legacy OTU dictionary and return a list of error handling results.

    If the OTU passes validation, an empty list is returned.

    :param fix: whether to attempt to fix the errors
    :param ncbi_client: the NCBI client
    :param otu: the otu data
    :return: a list of error handling results
    """
    try:
        LegacyOTU(**otu)
    except ValidationError as e:
        return handle_validation_error(e, fix, ncbi_client, otu)

    return None


def log_otu_validation_result(
    otu_name: str,
    otu_validation_result: OTUValidationResult,
    no_ok: bool,
) -> None:
    if not otu_validation_result.handler_results and no_ok:
        return

    console.rule(f"[bold]{otu_name}[/bold]", align="left", style="grey")
    console.line()

    if otu_validation_result.handler_results:
        table = Table.grid(padding=(0, 1, 0, 1))

        table.add_column()
        table.add_column()
        table.leading = 1

        for result in otu_validation_result.handler_results:
            label = "[green]FIXED[/green]" if result.fixed else "[red]ERROR[/red]"
            table.add_row(f"‣ [bold][{label}][/bold]", result.message)

        console.print(
            Padding(
                table,
                (0, 0, 0, 2),
            ),
        )
    else:
        console.print(
            Padding("‣ [green][bold]OK[/bold][/green]", (0, 0, 0, 2)),
        )

    console.line()


def validate_legacy_repo(fix: bool, limit: int, no_ok: bool, path: Path) -> None:
    """Validate a legacy reference repository."""
    ncbi = NCBIClient(ignore_cache=False)

    src_path = path / "src"

    with_errors_count = 0

    check_unique_accessions(path)
    check_unique_ids(path)
    check_unique_otu_abbreviations_and_names(path)

    for alpha_dir_path in sorted(src_path.iterdir()):
        if alpha_dir_path.is_file():
            continue

        for otu_dir_path in sorted(alpha_dir_path.iterdir()):
            if otu_dir_path.is_file():
                continue
            otu = build_legacy_otu(otu_dir_path)

            result = validate_legacy_otu(
                fix,
                ncbi,
                otu,
            )

            if result is None:
                continue

            if fix:
                replace_otu(
                    otu_dir_path,
                    result.repaired_otu,
                )

            log_otu_validation_result(
                result.repaired_otu["name"],
                result,
                no_ok,
            )

            with_errors_count += 1

            if with_errors_count >= limit:
                return

    logger.debug("Legacy repo passed validation.")
