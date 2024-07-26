"""User configuration path."""

from pathlib import Path

import platformdirs

user_cache_directory_path = Path(platformdirs.user_cache_dir("ref-builder"))
"""The path to a ref-builder specific cache directory in the user's system cache
path.

For Linux users, it would be ``$XDG_CACHE_HOME/ref_builder``.

"""
