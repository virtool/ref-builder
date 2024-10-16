from typing import Annotated, Union

from pydantic import Field, TypeAdapter

from ref_builder.events.repo import CreateRepo
from ref_builder.events.otu import (
    CreateOTU,
    CreatePlan,
    SetReprIsolate,
    ExcludeAccession,
)
from ref_builder.events.isolate import CreateIsolate, LinkSequence, DeleteIsolate
from ref_builder.events.sequence import CreateSequence, DeleteSequence


LoadableEvent = Annotated[
    Union[
        CreateRepo,
        CreateOTU,
        CreateIsolate,
        CreateSequence,
        DeleteIsolate,
        DeleteSequence,
        LinkSequence,
        CreatePlan,
        SetReprIsolate,
        ExcludeAccession,
    ],
    Field(discriminator="type"),
]


event_adapter = TypeAdapter(LoadableEvent)
