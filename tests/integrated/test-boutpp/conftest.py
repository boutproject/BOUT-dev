import pytest
import os
import sys
import subprocess


@pytest.fixture(autouse=True, scope="function")
def unique_xdist_group(request):
    # Unique group per test function (nodeid is unique, e.g., test_file.py::test_func)
    group_name = (
        f"boutpp_isolated_{request.node.nodeid.replace('/', '_').replace('::', '_')}"
    )
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
        sys.executable,
        "-m",
        "pytest",
        nodeid,
        "--rootdir",
        root_dir,
        "-p",
        "no:xdist",
        "-o",
        "cache_dir=/tmp/pytest-cache",  # Redirect cache to a writable place
        "-c",
        str(request.config.inifile or root_dir),  # Point to actual config if it exists
    ]

    result = subprocess.run(cmd, env=env, cwd=root_dir, capture_output=True, text=True)

    if result.returncode != 0:
        pytest.fail(
            f"Isolated test failed with exit code {result.returncode}\n"
            f"--- STDERR ---\n{result.stderr}\n"
            f"--- STDOUT ---\n{result.stdout}",
            pytrace=False,
        )

    # Return True so the parent can exit the test cleanly without skipping.
    return lambda: True


@pytest.fixture(autouse=True)
def sanitize_openmpi_env():
    """
    OpenMPI leaks state via environment variables (PMIX_*, OPAL_*).
    This fixture removes them so subprocesses called via `mpirun`
    don't crash thinking they are already initialized.
    """
    # Find all problematic OpenMPI state variables
    bad_prefixes = ("PMIX_", "OPAL_")
    bad_exact = (
        "OMPI_COMM_WORLD_SIZE",
        "OMPI_COMM_WORLD_RANK",
        "OMPI_COMM_WORLD_LOCAL_RANK",
    )

    mpi_vars_to_remove = [
        k for k in os.environ if k.startswith(bad_prefixes) or k in bad_exact
    ]

    # Save original variables
    saved_env = {k: os.environ[k] for k in mpi_vars_to_remove}

    # Delete them from the current environment
    for k in mpi_vars_to_remove:
        del os.environ[k]

    yield  # Run the test

    # Restore the environment after the test finishes
    os.environ.update(saved_env)
