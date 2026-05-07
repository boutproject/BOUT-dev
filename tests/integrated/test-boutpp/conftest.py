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
    Spawns a fresh, isolated process for tests.
    This prevents C++ singleton/MPI conflicts in xdist workers.
    Returns a function that evaluates to True in the parent (to abort test execution)
    and False in the child (to run the actual test).
    """
    if os.environ.get("BOUT_ISOLATED_RUN") == "1":
        return lambda: False

    # In the parent process, set up and run the child
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

    result = subprocess.run(cmd, env=env, cwd=root_dir, capture_output=True, text=True)

    if result.returncode != 0:
        pytest.fail(
            f"Isolated test failed with exit code {result.returncode}\n"
            f"--- STDERR ---\n{result.stderr}\n"
            f"--- STDOUT ---\n{result.stdout}",
            pytrace=False
        )

    # Return True so the parent can exit the test cleanly without skipping.
    return lambda: True
