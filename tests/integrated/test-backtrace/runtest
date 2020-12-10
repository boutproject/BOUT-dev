#!/usr/bin/env python3

# Test enabling/disabling exception backtrace from environment variable

# requires all_tests

from boututils.run_wrapper import build_and_log, shell
import os

build_and_log("backtrace environment variable test")

try:
    del os.environ["BOUT_SHOW_BACKTRACE"]
except KeyError:
    pass

success = True
errors = []

_, output = shell("./boutexcept", pipe=True)

if "troublemaker" in output:
    success = False
    print("Fail: detected offending function name in output when not expected")

_, output = shell("BOUT_SHOW_BACKTRACE=yes ./boutexcept", pipe=True)

if "troublemaker" not in output:
    success = False
    print("Fail: did not detect offending function name in output when expected")

if success:
    print("=> All BoutException backtrace tests passed")
    exit(0)

print("=> Some failed tests")
exit(1)
