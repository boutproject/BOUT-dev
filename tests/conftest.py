import os
import pytest
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def chdir(path: Path):
    """Temporarily change the current working directory."""
    old_cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_cwd)

@pytest.fixture
def test_dir(request) -> Path:
    return Path(request.fspath).parent
