from pydantic import UUID4

from ref_builder.events.base import EventData, Event, RepoQuery


class CreateRepoData(EventData):
    """The data for a :class:`CreateRepo` event."""

    id: UUID4
    data_type: str
    name: str
    organism: str


class CreateRepo(Event):
    """An event that creates a new repository.

    This event is always the first event in a repository's event log.
    """

    data: CreateRepoData
    query: RepoQuery
