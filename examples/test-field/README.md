test-field
==========

Test that arithmetic operations on fields give the correct results.

This test initialises a Field2D and a Field3D using BoutReals, then performs
various arithmetic operations on them, checking that different indices are
equal to the correct result with a tolerance of 1e-10.

The check is performed using the ASSERT0 macro (defined in assert.hxx).

The following operators are tested:

- Field2D:
    - `Field2D::operator+=(const BoutReal)`
    - `Field2D::operator-=(const BoutReal)`
    - `Field2D::operator*=(const BoutReal)`
    - `Field2D::operator/=(const BoutReal)`

    - `Field2D::operator+=(const Field2D&)`
    - `Field2D::operator*=(const Field2D&)`

- Field3D:
    - `Field3D::operator+=(const BoutReal)`
    - `Field3D::operator-=(const BoutReal)`
    - `Field3D::operator*=(const BoutReal)`
    - `Field3D::operator/=(const BoutReal)`

    - `Field3D::operator-=(const Field2D&)`
    - `Field3D::operator*=(const Field2D&)`

    - `Field3D::operator+=(const Field3D&)`
    - `Field3D::operator*=(const Field3D&)`
    - `Field3D::operator/=(const Field3D&)`

Additionally, the non-member function `pow(const Field3D &lhs, const Field2D
&rhs)` is tested.
