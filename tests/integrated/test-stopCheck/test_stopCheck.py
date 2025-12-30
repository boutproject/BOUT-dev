import os
import subprocess
import time
from contextlib import chdir


def test_runtest(test_dir):
    # BOUT_TOP from script dir (CWD-independent)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    bout_root = os.path.abspath(os.path.join(script_dir, ".."))
    os.environ["BOUT_TOP"] = bout_root

    # Add pylib (dependencies like boutdata) to PYTHONPATH
    pylib_path = os.path.join(bout_root, "tools", "pylib")
    if os.path.exists(pylib_path):
        os.environ["PYTHONPATH"] = f"{pylib_path}:{os.environ.get('PYTHONPATH', '')}"

    # Pre-build for dependencies like grid.fci.nc

    with chdir(test_dir):
        build_result = subprocess.run("make", shell=True, cwd=test_dir, capture_output=True, text=True, timeout=60)
        if build_result.returncode != 0:
            print(f"Build stderr in {test_dir}: {build_result.stderr}")  # Non-fatal; some tests no-op

        # MPI oversubscribe for communications test
        cmd = './runtest'
        if "communications" in str(test_dir):
            os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "1"  # Allows 18 procs
        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
        elapsed = time.time() - start
        print(f"Output: {result.stdout}")
        if result.returncode != 0:
            print(f"Stderr: {result.stderr}")
        assert result.returncode == 0, f"Failed after {elapsed:.3f}s in {test_dir}"
