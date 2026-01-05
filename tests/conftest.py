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
    Provides a unique, writable copy of a pre-existing input directory
    and changes cwd to it.

    The original directory name is taken from:
      - or @pytest.mark.input_dir("name") marker
      - default: "data"
    """

    marker = request.node.get_closest_marker("input_dir")
    original_name = marker.kwargs.get('name', 'data') if marker else "data"

    test_file_dir = Path(request.fspath).parent
    original_dir = test_file_dir / original_name

    if not original_dir.is_dir():
        pytest.fail(f"Expected input directory '{original_dir}' not found")

    # Create unique writable copy
    run_dir_path = tmp_path / original_name
    shutil.copytree(original_dir, run_dir_path, dirs_exist_ok=True)

    monkeypatch.chdir(tmp_path)

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
