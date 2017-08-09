Test FieldGroup class
=====================

The FieldGroup class stores a collection of fields, and is used in communications. 
Variables stored in FieldGroup should be derived from FieldData.

This test consists of two parts

1. Compile a short code (`test_fail.cxx`) which should fail to compile.
   This code tries to construct a FieldGroup by passing an int.

2. A test which adds some fields to a FieldGroup then checks that
   the correct number of fields are present.


