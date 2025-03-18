from ref_builder.console import console, print_otu_event_log
from ref_builder.repo import Repo


def test_otu_event_log(scratch_repo: Repo):
    otu = next(scratch_repo.iter_otus())

    otu_events = list(scratch_repo.iter_otu_events(otu.id))

    with console.capture() as capture:
        print_otu_event_log(otu_events)

