import pytest
from pathlib import Path


@pytest.fixture(autouse=True)
def cwd_to_test_file_dir(monkeypatch, test_dir):
    """Automatically change CWD to the directory of the current test file."""
    monkeypatch.chdir(test_dir)

@pytest.fixture
def test_dir(request) -> Path:
    return Path(request.fspath).parent
