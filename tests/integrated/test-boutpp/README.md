Testing the boutpp interface
=============================

Some tests for the boutpp interface:


collect
-------
Basic test wether the fromCollect is working

collect-staggered
-----------------
Test on whether we are reading the staggering metadata correctly, and
setting it correctly. It further verifies that interp_to on the so
loaded fields in interpolating correctly - so that loading staggered
fields and using them for other locations works as intended.

legacy-model
------------
A test that the old interface is not crashing

mms-ddz
-------
A simple convergence test for DDZ and C2.
Mostly there to make sure DDZ returns a derivative

simple-model
-----------
Testing whether we can derive from the PhysicsModel without crashing.
