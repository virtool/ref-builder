from collections.abc import Iterator

import rich.console
from rich.table import Table
from rich.text import Text

from ref_builder.models import OTUMinimal
from ref_builder.plan import SegmentRule
from ref_builder.resources import RepoOTU


def _render_taxonomy_id_link(taxid: int) -> str:
    return f"[link=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}]{taxid}[/link]"


def _render_nucleotide_link(accession: str) -> str:
    return f"[link=https://www.ncbi.nlm.nih.gov/nuccore/{accession}]{accession}[/link]"


def print_otu_as_json(otu: RepoOTU) -> None:
    """Print the OTU data to the console as JSON."""
    console.print(otu.model_dump_json())


def print_otu(otu: RepoOTU) -> None:
    """Print the details for an OTU to the console.

    :param otu: The OTU to print.

    """
    console.print(Text(otu.name, style="bold underline"))
    console.line()

    table = Table(
        box=None,
        show_header=False,
    )

    table.add_row("[bold]ACRONYM[/bold]", otu.acronym)
    table.add_row("[bold]ID[/bold]", str(otu.id))

    if otu.legacy_id:
        table.add_row("[bold]LEGACY ID[/bold]", otu.legacy_id)

    table.add_row("[bold]TAXID[/bold]", _render_taxonomy_id_link(otu.taxid))

    max_accession_length = max(
        len(str(sequence.accession))
        for isolate in otu.isolates
        for sequence in isolate.sequences
    )

    max_segment_name_length = max(
        len(str(segment.name)) for segment in otu.plan.segments
    )

    console.print(table)

    console.line()
    console.print("[bold]PLAN[/bold]")
    console.line()

    plan_table = Table(box=None)

    plan_table.add_column("NAME")
    plan_table.add_column("REQUIRED")
    plan_table.add_column("LENGTH")
    plan_table.add_column("TOLERANCE")
    plan_table.add_column("ID")

    for segment in otu.plan.segments:
        plan_table.add_row(
            str("Unnamed" if segment.name is None else segment.name),
            "[red]Yes[/red]"
            if segment.required == SegmentRule.REQUIRED
            else "[grey]No[/grey]",
            str(segment.length),
            str(segment.length_tolerance),
            str(segment.id),
        )

    console.print(plan_table)

    console.line()
    console.print("[bold]ISOLATES[/bold]")

    index_by_segment_id = {segment.id: i for i, segment in enumerate(otu.plan.segments)}

    for isolate in otu.isolates:
        console.line()
        console.print(
            Text(str(isolate.name))
        ) if isolate.name is not None else console.print(
            "[italic]Unnamed[/italic]",
        )
        console.line()

        isolate_table = Table(
            box=None,
        )

        isolate_table.add_column("ACCESSION", width=max_accession_length)
        isolate_table.add_column("LENGTH")
        isolate_table.add_column("SEGMENT", min_width=max_segment_name_length)
        isolate_table.add_column("DEFINITION")

        for sequence in sorted(
            isolate.sequences,
            key=lambda s: index_by_segment_id[s.segment],
        ):
            isolate_table.add_row(
                _render_nucleotide_link(str(sequence.accession)),
                str(len(sequence.sequence)),
                str(otu.plan.get_segment_by_id(sequence.segment).name or "Unnamed"),
                sequence.definition,
            )

        console.print(isolate_table)


def print_otu_list(otus: Iterator[OTUMinimal]) -> None:
    """Print a list of OTUs to the console.

    :param otus: The OTUs to print.

    """
    rows = [(otu.name, otu.acronym, str(otu.taxid), str(otu.id)) for otu in otus]

    if not rows:
        console.print("No OTUs found")
        return

    raw_headers = ("NAME", "ACRONYM", "TAXID", "ID")
    col_widths = [max(len(row[i]) for row in rows + [raw_headers]) for i in range(4)]

    # Define a row format based on column widths
    row_format = "  ".join(f"{{:<{w}}}" for w in col_widths)

    # Header without bold.
    header_text = row_format.format(*raw_headers)

    # Apply ANSI bold code and print.
    print(f"\033[1m{header_text}\033[0m")  # noqa: T201
    print("\n".join([row_format.format(*row) for row in rows]))  # noqa: T201


console = rich.console.Console()
