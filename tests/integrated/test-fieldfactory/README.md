test-fieldfactory
=================

Test creating Field2Ds and Field3Ds using a FieldFactory and compare the output
to existing benchmark files, using different numbers of processes.

This test is different to test-initial in that this test checks the ability of
calling the `FieldFactory::create2D` and `FieldFactory::create3D` methods
directly in a physics model, whereas test-initial is designed to check the
actual FieldGenerators.

This test is run using 1, 2 and 4 MPI processes.
