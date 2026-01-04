import pytest
import shutil
import boutpp
from pathlib import Path


@pytest.fixture(autouse=True)
def cwd_to_test_file_dir(monkeypatch, test_dir):
    """Automatically change CWD to the directory of the current test file."""
    monkeypatch.chdir(test_dir)

@pytest.fixture
def test_dir(request) -> Path:
    return Path(request.fspath).parent

@pytest.fixture
def copy_data_directory(request, tmp_path):
    """
    Provides a writable copy of a pre-existing directory (e.g. "data", "test", "input")
    in a unique temporary location.

    The name of the original directory is taken from request.node.input_dir_name
    (set via a marker on the test).
    """
    # Get the requested original directory name (default to "data" if not specified)
    original_name = getattr(request.node, "input_dir_name", "data")

    # Location of the original input directory (sibling to the test file)
    test_file_dir = Path(request.fspath).parent
    original_dir = test_file_dir / original_name

    if not original_dir.is_dir():
        pytest.fail(f"Expected input directory '{original_dir}' does not exist")

    # Unique writable copy for this test
    run_dir = tmp_path / original_name
    shutil.copytree(original_dir, run_dir, dirs_exist_ok=True)

    return run_dir

@pytest.fixture(autouse=True)
def finalize_boutpp(request):
    """Automatically finalize BOUT++ after each test that used it."""
    yield  # Run the test

    # Check if this test used boutpp init (heuristic: look for marker or direct check)
    if hasattr(request.node, "boutpp_initialized"):
        try:
            boutpp.finalise()
        except Exception:
            pass  # If already finalized or error, ignore
