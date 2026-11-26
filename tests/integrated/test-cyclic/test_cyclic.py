#!/usr/bin/env python3

import pytest
from boututils.run_wrapper import launch_safe


@pytest.mark.parametrize("flags", ["", "nsys=2", "nsys=5 periodic", "nsys=7 n=10"])
@pytest.mark.parametrize("nproc", [1, 2, 4])
def test_cyclic(flags, nproc):
    launch_safe(f"./test_cyclic {flags}", nproc=nproc, mthread=1)
