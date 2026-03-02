#!/usr/bin/env python3

# requires: all_tests

from boututils.run_wrapper import shell, shell_safe

from boutdata import collect
import numpy as np


def test_collect():


    shell_safe("./test-collect")

    # Try collecting data using incorrect case
    # This should be corrected automatically
    a = collect("A", path="data")

    assert np.allclose(a, 1.23), f"Wrong value => Failed: {a}"
