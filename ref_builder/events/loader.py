from datetime import datetime
from enum import StrEnum
from typing import Annotated, Union

from pydantic import BaseModel, Field, TypeAdapter

from ref_builder.events.repo import CreateRepo
from ref_builder.events.otu import CreateOTU, CreatePlan, SetReprIsolate, ExcludeAccession
from ref_builder.events.isolate import CreateIsolate, LinkSequence, DeleteIsolate
from ref_builder.events.sequence import CreateSequence, DeleteSequence


class EventType(StrEnum):
    CREATE_ISOLATE = "CreateIsolate"
    CREATE_OTU = "CreateOTU"
    CREATE_PLAN = "CreatePlan"
    CREATE_REPO = "CreateRepo"
    CREATE_SEQUENCE = "CreateSequence"
    DELETE_ISOLATE = "DeleteIsolate"
    DELETE_SEQUENCE = "DeleteSequence"
    EXCLUDE_ACCESSION = "ExcludeAccession"
    LINK_SEQUENCE = "LinkSequence"
    SET_REPR_ISOLATE = "SetReprIsolate"


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
    Field(discriminator='type'),
]


event_adapter = TypeAdapter(LoadableEvent)