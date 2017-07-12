Restart test
============

This checks that code restarting works, by running three cases:

* A simulation with 10 outputs

* A simulation with 5 outputs, restart and appending another 5 outputs.
  The result should be the same as the first run.

* Run 5 outputs, then restart without appending for another 5 outputs.
  This should contain the initial value as first time

