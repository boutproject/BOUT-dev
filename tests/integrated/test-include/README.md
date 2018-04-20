test-include
============

This test ensures that all header files from bout can be included on
there own. It compiles a test file with one header, and if that fails,
then the headerfile does not include all files it needs.

Checking whether the library builds as whole is not sufficient, as
ether the file might never be included directly, or some other header
file might always be included first, which fullfills the requirements,
thus there is this extra test.
