#!/usr/bin/env python3

# Test enabling/disabling exception backtrace from environment variable

from boututils.run_wrapper import shell
import os


def test_backtrace():
    try:
        del os.environ["BOUT_SHOW_BACKTRACE"]
    except KeyError:
        pass

    _, output = shell(["BOUT_SHOW_BACKTRACE=0 ./boutexcept"], pipe=True)

    assert "troublemaker" not in output, (
        f"Fail: detected offending function name in output when not expected\n{output}"
    )

    _, output = shell(["./boutexcept"], pipe=True)

    assert "troublemaker" in output, (
        f"Fail: did not detect offending function name in output when expected\n{output}"
    )
