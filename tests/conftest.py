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
