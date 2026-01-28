import pytest
import shutil
import boutpp
from pathlib import Path


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "input_dir(name): specify the pre-existing input directory name for this test"
    )
    config.addinivalue_line(
        "markers", "serial: mark that the test should not be run concurrently with other. "
    )

@pytest.fixture
def test_dir(request) -> Path:
    return Path(request.fspath).parent

@pytest.fixture(scope="function", autouse=True)
def copy_and_cwd_to_unique_tmp_dir(request, tmp_path_factory, monkeypatch):
    """
    For each test function, create a unique temporary copy of the test directory
    and change cwd to it.
    """

    test_file_dir = Path(request.fspath).parent

    if not test_file_dir.is_dir():
        pytest.fail(f"Expected test directory '{test_file_dir}' not found")

    # Create a unique temp dir for this test
    run_dir = tmp_path_factory.mktemp(test_file_dir.name)

    # Copy the original test directory into it
    shutil.copytree(test_file_dir, run_dir, dirs_exist_ok=True)

    # Change working directory to the copy
    monkeypatch.chdir(run_dir)
