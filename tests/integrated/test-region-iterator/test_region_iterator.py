#!/usr/bin/env python3

#
# Run the test, check it completed successfully
#

from boututils.run_wrapper import launch_safe, shell

flags = [""]
cmd = "./test_region_iterator"
code = 0  # Return code
pipe = True


def test_region_iterator():

    for nproc in [1, 2]:  # Number of mpi procs
        for mthread in [1]:  # Number of omp threads (not yet supported)
            print(f"\t{nproc} processors and {mthread} threads")

            for f in flags:
                # Run the case
                shell(["rm -rf data/BOUT.dmp.*"])
                _s, out = launch_safe(cmd + " " + f, nproc=nproc, pipe=pipe)
                if pipe:
                    with open(f"run.log.{nproc}.{mthread}", "w") as f:
                        f.write(out)

                # If we've got here we know that cmd launched by launch_safe passed
                # as otherwise it raises.
                print("Passed")

    assert code == 0, " => Some failed tests"
