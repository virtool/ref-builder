import shutil
from pathlib import Path

import orjson
import pytest
from pytest_mock import MockerFixture

from ref_builder.legacy.utils import build_legacy_otu
from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.client import NCBIClient
from ref_builder.repo import Repo
from ref_builder.utils import DataType
from tests.src_test_gen import generate_src_test_files


@pytest.fixture()
def uncached_ncbi_client(scratch_ncbi_cache: NCBICache) -> NCBIClient:
    """An NCBI client that ignores the cache."""
    scratch_ncbi_cache.clear()
    return NCBIClient(ignore_cache=True)


@pytest.fixture()
def files_path():
    return Path(__file__).parent / "files"


@pytest.fixture()
def empty_repo(tmp_path: Path) -> Repo:
    """An empty reference repository."""
    return Repo.new(
        DataType.GENOME,
        "Generic Viruses",
        tmp_path / "test_repo",
        "virus",
    )


@pytest.fixture()
def legacy_otu(legacy_repo_path: Path) -> dict:
    """A legacy OTU."""
    return build_legacy_otu(
        legacy_repo_path / "src" / "a" / "abaca_bunchy_top_virus",
    )


@pytest.fixture()
def legacy_repo_path(
    files_path: Path,
    tmp_path: Path,
) -> Path:
    """The path to a scratch legacy reference."""
    path = tmp_path / "legacy_repo"

    shutil.copytree(files_path / "src_v1", path / "src", dirs_exist_ok=True)

    return path


@pytest.fixture()
def precached_repo(
    mocker: MockerFixture,
    scratch_user_cache_path: Path,
    tmp_path: Path,
) -> Repo:
    """A reference repository with a preloaded NCBI cache."""
    mocker.patch(
        "ref_builder.paths.user_cache_directory_path",
        return_value=scratch_user_cache_path,
    )

    return Repo.new(DataType.GENOME, "Empty", tmp_path / "precached_repo", "virus")


@pytest.fixture()
def scratch_ncbi_cache(mocker: MockerFixture, scratch_user_cache_path: Path):
    """A scratch NCBI cache with preloaded data."""
    mocker.patch(
        "ref_builder.ncbi.cache.user_cache_directory_path",
        scratch_user_cache_path,
    )

    return NCBICache()


@pytest.fixture()
def scratch_ncbi_client(
    mocker: MockerFixture,
    scratch_user_cache_path: Path,
) -> NCBIClient:
    """A scratch NCBI client with a preloaded cache."""
    mocker.patch(
        "ref_builder.ncbi.cache.user_cache_directory_path",
        scratch_user_cache_path,
    )

    return NCBIClient(ignore_cache=False)


@pytest.fixture()
def scratch_path(scratch_repo: Repo) -> Path:
    """The path to a scratch reference repository."""
    return scratch_repo.path


@pytest.fixture()
def scratch_repo(files_path: Path, tmp_path: Path, scratch_event_store_data: dict[str, dict]) -> Repo:
    """A prepared scratch repository."""
    path = tmp_path / "scratch_repo"

    src_path = path / "src"
    src_path.mkdir(parents=True)

    for filename in scratch_event_store_data:
        with open(src_path / filename, "wb") as f:
            f.write(orjson.dumps(scratch_event_store_data[filename]))

    (path / ".cache").mkdir(parents=True)

    return Repo(path)


@pytest.fixture()
def scratch_user_cache_path(files_path: Path, tmp_path: Path) -> Path:
    """A path to a user cache that contains preloaded data.

    On Linux, a user cache path would be found at ``~/.cache/ref-builder/ncbi``.

    """
    path = tmp_path / "user_cache"
    path.mkdir()

    shutil.copytree(files_path / "cache_test", path / "ncbi")

    return path


@pytest.fixture()
def scratch_repo_contents_path(files_path):
    return files_path / "src_test_contents.json"


@pytest.fixture
def scratch_event_store_data(pytestconfig, tmp_path, scratch_repo_contents_path):
    scratch_src = pytestconfig.cache.get("scratch_src", {})
    if not scratch_src:
        scratch_src = generate_src_test_files(tmp_path, scratch_repo_contents_path)
        if scratch_src:
            pytestconfig.cache.set("scratch_src", scratch_src)
        else:
            raise (ValueError("scratch_src did not generate properly"))

    return scratch_src