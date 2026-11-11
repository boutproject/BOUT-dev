#!/usr/bin/env python3

# Run the test, compare results against the benchmark

import sys

from boutdata.collect import collect
from boututils.run_wrapper import launch

nproc = 1


def test_stopCheck_file():

    check_values = [True, False]
    expected_steps = [1, 11]
    data_dirs = ["data", "dataSecond"]
    stop_files = ["BOUT.stop", "otherStop.check"]

    print("Running stopCheck test")
    success = True

    for data_dir, stop_file in zip(data_dirs, stop_files):
        for check, expected in zip(check_values, expected_steps):
            executable = "./test_stopCheck"
            command = f"{executable} -d {data_dir} stopCheck={check} stopCheckName={stop_file} replace"

            _s, out = launch(
                command,
                nproc=nproc,
                pipe=True,
            )
            with open(f"run.log.{data_dir}.{check}", "w") as f:
                f.write(out)

            result = collect("t_array", path=data_dir, info=False)

            if result.shape[0] != expected:
                print("Fail, wrong shape")
                print(f"\tOption is {check}/{data_dir}/{stop_file}")
                print(f"\tshape is {result.shape[0]}")
                print(f"\texpecting {expected}")
                sys.exit(1)
                success = False
                continue

    assert success, " => Some failed tests"
