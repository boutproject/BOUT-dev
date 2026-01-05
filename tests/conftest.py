import pytest
import shutil
import boutpp
from pathlib import Path


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "input_dir(name): specify the pre-existing input directory name for this test"
    )

@pytest.fixture
def test_dir(request) -> Path:
    return Path(request.fspath).parent

@pytest.fixture(autouse=True)
def copy_and_cwd_to_tmp_dir(request, tmp_path, monkeypatch):
    """
    Provides a unique, writable copy of the test directory
    and changes cwd to it.
    """

    test_file_dir = Path(request.fspath).parent

    if not test_file_dir.is_dir():
        pytest.fail(f"Expected input directory '{test_file_dir}' not found")

    # Create unique writable copy
    run_dir_path = tmp_path / test_file_dir.name
    shutil.copytree(test_file_dir, run_dir_path, dirs_exist_ok=True)

    monkeypatch.chdir(run_dir_path)
