import pytest
import os
import sys
import subprocess


@pytest.fixture(autouse=True, scope="function")
def unique_xdist_group(request):
    # Unique group per test function (nodeid is unique, e.g., test_file.py::test_func)
    group_name = f"boutpp_isolated_{request.node.nodeid.replace('/', '_').replace('::', '_')}"
    request.node.add_marker(pytest.mark.xdist_group(name=group_name))

@pytest.fixture(scope="function")
def run_isolated(request):
    """
    Fixture to run a test in a fresh, isolated subprocess.
    This prevents C++ singleton/MPI conflicts in xdist workers.
    """
    # If this variable is set, we are already inside the isolated subprocess
    if os.environ.get("BOUT_ISOLATED_RUN") == "1":
        return

    # Get the unique ID of the current test (e.g., path/to/test.py::test_func)
    nodeid = request.node.nodeid

    # Use the absolute path to the test file to avoid directory resolution issues
    # fspath is the local path to the test file
    root_dir = str(request.config.rootpath)

    env = os.environ.copy()
    env["BOUT_ISOLATED_RUN"] = "1"

    # Disable xdist (-p no:xdist) in the subprocess to ensure a clean single-process run
    # Use -o to disable the cache entirely in the subprocess
    # Remove -c /dev/null so it stays in the project context
    cmd = [
        sys.executable, "-m", "pytest",
        nodeid,
        "--rootdir", root_dir,
        "-p", "no:xdist",
        "-o", "cache_dir=/tmp/pytest-cache", # Redirect cache to a writable place
        "-c", str(request.config.inifile or root_dir) # Point to actual config if it exists
    ]

    # Run from the same directory as the parent
    result = subprocess.run(
        cmd,
        env=env,
        cwd=root_dir,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        pytest.fail(
            f"Isolated test failed with exit code {result.returncode}\n"
            f"--- STDERR ---\n{result.stderr}\n"
            f"--- STDOUT ---\n{result.stdout}",
            pytrace=False
        )

    # If the subprocess succeeded, skip the execution of the test body
    # in the parent (original) xdist worker.
    pytest.skip("Test successfully completed in isolated subprocess")
