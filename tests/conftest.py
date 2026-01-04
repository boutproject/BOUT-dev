import pytest
import boutpp
from pathlib import Path


@pytest.fixture(autouse=True)
def cwd_to_test_file_dir(monkeypatch, test_dir):
    """Automatically change CWD to the directory of the current test file."""
    monkeypatch.chdir(test_dir)

@pytest.fixture
def test_dir(request) -> Path:
    return Path(request.fspath).parent

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
