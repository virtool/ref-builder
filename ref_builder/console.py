from collections.abc import Iterable

import rich.console
from rich.table import Table

from ref_builder.resources import RepoOTU


def _render_taxonomy_id_link(taxid: int) -> str:
    return f"[link=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}]{taxid}[/link]"


def print_otu(otu: RepoOTU) -> None:
    """Print the details for an OTU to the console.

    :param otu: The OTU to print.

    """
    console.print(f"[bold][underline]{otu.name}[/bold][/underline]")
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
        len(sequence.accession)
        for isolate in otu.isolates
        for sequence in isolate.sequences
    )

    max_segment_name_length = max(len(segment.name) for segment in otu.schema.segments)

    console.print(table)

    console.line()
    console.print("[bold]SCHEMA[/bold]")
    console.line()

    schema_table = Table(box=None)

    schema_table.add_column("SEGMENT")
    schema_table.add_column("REQUIRED")
    schema_table.add_column("LENGTH")

    for segment in otu.schema.segments:
        schema_table.add_row(
            segment.name,
            "[red]Yes[/red]" if segment.required else "[grey]No[/grey]",
            str(segment.length),
        )

    console.print(schema_table)

    console.line()
    console.print("[bold]ISOLATES[/bold]")

    for isolate in otu.isolates:
        console.line()
        console.print(str(isolate.name))
        console.line()

        isolate_table = Table(
            box=None,
        )

        isolate_table.add_column("ACCESSION", width=max_accession_length)
        isolate_table.add_column("SEGMENT", min_width=max_segment_name_length)
        isolate_table.add_column("DEFINITION")

        for sequence in sorted(isolate.sequences, key=lambda s: s.accession):
            isolate_table.add_row(
                sequence.accession,
                sequence.segment,
                sequence.definition,
            )

        console.print(isolate_table)


def print_otu_list(otus: Iterable[RepoOTU]) -> None:
    """Print a list of OTUs to the console.

    :param otus: The OTUs to print.

    """
    table = Table(
        box=None,
    )

    table.add_column("NAME")
    table.add_column("ACRONYM")
    table.add_column("TAXID")
    table.add_column("ID")

    for otu in otus:
        table.add_row(otu.name, otu.acronym, str(otu.taxid), str(otu.id))

    console.print(table)


console = rich.console.Console()
