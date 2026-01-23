import os
import subprocess
import time


def test_runtest(test_dir):
    cmd = './runtest'
    # MPI oversubscribe for communications test
    os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "1"  # Allows 18 procs
    start = time.time()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
    elapsed = time.time() - start
    print(f"Output: {result.stdout}")
    if result.returncode != 0:
        print(f"Stderr: {result.stderr}")
    assert result.returncode == 0, f"Failed after {elapsed:.3f}s in {test_dir}"
