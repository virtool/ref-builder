from uuid import UUID

from ref_builder.resources import RepoOTU
from ref_builder.utils import IsolateName, IsolateNameType


def get_isolate_id_from_key(otu: RepoOTU, isolate_key: str):
    try:
        return UUID(isolate_key)
    except ValueError:
        parts = isolate_key.split(" ")
        try:
            isolate_name = IsolateName(IsolateNameType(parts[0].lower()), parts[1])
            return otu.get_isolate_id_by_name(isolate_name)
        except ValueError:
            raise ValueError(f'Error: "{isolate_key}" is not a valid isolate name.')
