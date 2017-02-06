Notes
=====

Compile options
---------------

Compiling with ``-DCHECK`` enables a lot of checks of operations
performed by the field objects. This is very useful for debugging a
code, and can be omitted once bugs have been removed.

For (sometimes) more useful error messages, there is the ``-DTRACK``
option. This keeps track of the names of variables and includes these in
error messages.

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

