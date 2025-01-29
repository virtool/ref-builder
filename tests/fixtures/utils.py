from uuid import UUID

from syrupy.matchers import path_type

uuid_matcher = path_type({".*id": (UUID,)}, regex=True)
"""A Syrupy matcher for UUIDs.

Matches any key ending with "id" to a UUID object.
"""
