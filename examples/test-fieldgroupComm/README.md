test-fieldgroupComm
===================

Test communicating FieldGroups for different number of processes, checking the
results against a "correct" answer.

Three identical Field3Ds are created and added in different combinations to
three separate communicators. One communicator is used "correctly" and is
defined as giving the correct answer; the second contains two copies of the same
field, and the third is communicated twice in a row. `Grad_par` is then called
on the fields.

The results of the second and third fields are compared against the first with a
tolerance of 1e-10.
