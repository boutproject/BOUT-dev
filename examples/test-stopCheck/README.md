Stop check test case
====================

This example tests the `stopCheck` functionality by checking that the number of output steps matches the expectations for a run with `stopCheck=false` and one with `stopCheck=true`.

Note the data*/BOUT.stop file is expected to exist and is not created by the test runner.

This tests both the behaviour of `checkStop` with default options as well as the behaviour with non-default data_dir and stop file name.