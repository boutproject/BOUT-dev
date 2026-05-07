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

    # Copy current environment and set the isolation flag
    env = os.environ.copy()
    env["BOUT_ISOLATED_RUN"] = "1"

    # Construct the command to run ONLY this specific test
    # We disable xdist (-p no:xdist) in the subprocess to ensure a clean single-process run
    cmd = [
        sys.executable, "-m", "pytest",
        "-p", "no:xdist",
        "-c", "/dev/null", # Ignore config to prevent recursive loading if necessary
        nodeid
    ]

    result = subprocess.run(cmd, env=env, capture_output=True, text=True)

    if result.returncode != 0:
        # If the isolated run failed, fail the parent test with the output
        pytest.fail(
            f"Isolated test failed with exit code {result.returncode}\n"
            f"--- STDERR ---\n{result.stderr}\n"
            f"--- STDOUT ---\n{result.stdout}",
            pytrace=False
        )

    # If the subprocess succeeded, we skip the execution of the test body
    # in the parent (original) xdist worker.
    pytest.skip("Test successfully completed in isolated subprocess")
