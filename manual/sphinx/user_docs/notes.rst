Notes
=====

Configure options
-----------------

Configure with ``--enable-checks=3`` enables a lot of checks of
operations performed by the field objects. This is very useful for
debugging a code, and can be omitted once bugs have been removed.
``--enable=checks=2`` enables less checking, especially the
computationally rather expensive ones, while ``--enable-checks=0``
disables most checks.

To get most checking, both from BOUT++ and from the compiler
``--enable-debug`` can be used. That enables checks of level 3, as
well as debug flags, e.g. ``-g`` for gcc.

For (sometimes) more useful error messages, there is the
``--enable-track`` option. This keeps track of the names of variables
and includes these in error messages.

To enable optimization, configure with ``--enable-optimize=3``.
This will try to set appropriate flags, but may not set the best ones.
This should work well for gcc. Similar to checks, different levels can
be specified, where 3 is high, and 0 means disabling all
optimization. ``--enable-optimize=fast`` will set the ``-Ofast`` flag
for gcc which enables optimization that are not standard conform, so
proceed at own risk.

Compile options
---------------

It is possible to change flags for BOUT++ after running configure, by
editing the make.config file. Note that this is not recommended, as
e.g. pvode will not be build with these flags.

Adaptive grids
--------------

Two types of adaptive grids can be used in BOUT++: Moving meshes, and
changing resolution.

Moving meshes
~~~~~~~~~~~~~

During either the initialisation, or the simulation itself, the metric
tensors can be modified. This could be used to make the coordinate
system time-dependent. Since currently the metric tensors are 2D fields,
this would only allow axisymmetric motion. Changing the tensors to be 3D
objects is however possible with fairly small modification to the code.

Whenever one of the metrics :math:`g^{ij}` are changed, a call to
``geometry()`` must be made.

Changing resolution
~~~~~~~~~~~~~~~~~~~

Since all 2D and 3D fields/vectors are located internally in global
lists, the resolution of the grid can be changed when required by
interpolation. **This requires a new, more efficient implementation of
the Fields classes**.

