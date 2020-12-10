blob2d using `LaplaceXZ`
========================

This is very similar to the 2D drift-reduced [blob2d
model](../examples/README.md), except that the perpendicular Laplacian inversion
solver for `phi` uses [`LaplaceXZ`][laplacexz] instead of
[`Laplace`][laplace]. See the linked documentation for details on the
differences


[laplacexz]: https://bout-dev.readthedocs.io/en/latest/user_docs/laplacian.html#laplacexz
[laplace]: https://bout-dev.readthedocs.io/en/latest/user_docs/laplacian.html#laplacian-inversion
