.. _sec-iterating:

Iterating over fields
=====================

As of BOUT++ 4.0.0, we now have the ability to use C++ range-based
for-loops. This means that it is possible to iterate over a whole field
using a single loop:

::

    Field3D f(0.0);
    for (auto i : f) {
       f[i] = a[i] + b[i];
    }

This replaces the C-style triple-nested loop:

::

   Field3D f(0.0);
   for (int i = mesh->xstart; i < mesh->xend; ++i) {
     for (int j = mesh->ystart; j < mesh->yend; ++j) {
       for (int k = 0; k < mesh->LocalNz; ++k) {
         f[i,j,k] = a[i,j,k] + b[i,j,k]
       }
     }
   }

The iterator provides access to the x, y, z indices:

::

    Field3D f(0.0);
    for (auto i : f) {
       f[i] = i.x + i.y + i.z;
    }

It is also possible to specify regions to iterate over using this
syntax:

::

    Field3D f(0.0);
    for (auto i : f.region(RGN_NOBNDRY)) {
       f[i] = 1.0;
    }

Available regions are:

-  ``RGN_ALL``, which is the whole mesh;

-  ``RGN_NOBNDRY``, which skips all boundaries;

-  ``RGN_NOX``, which skips the x boundaries

-  ``RGN_NOY``, which skips the y boundaries
