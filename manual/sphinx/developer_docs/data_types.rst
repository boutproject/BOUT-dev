Arrays, scalar and vector field types
=====================================

The classes outlines in red in :numref:`fig-layout2` are data types
currently implemented in BOUT++.

``FieldData``
-------------

All BOUT++ data types implement a standard interface for accessing their
data, which is then used in communication and file I/O code. This
interface is in ``src/field/field_data.hxx``. The mandatory (pure
virtual) functions are:

::

    bool isReal(); // Returns true if field consists of real values
    bool is3D() const;   // True if variable is 3D
      
    int byteSize() const; // Number of bytes for a single point
    int realSize() const; // Number of reals (not implemented if not real)

To support file I/O there are also some additional functions which may
be implemented. A code can check if they are implemented by calling
``ioSupport``. If one of them is implemented then they all should be.

::

    bool  ioSupport();  // Return true if these functions are implemented
    const string getSuffix(int component) const; // For vectors e.g. "_x"
    void* getMark() const; // Store current settings (e.g. co/contra-variant)
    void  setMark(void *setting); // Return to the stored settings
    BoutReal* getData(int component); 
    void  zeroComponent(int component); // Set a component to zero

For twist-shift conditions, the optional function ``shiftZ`` is called
in the communication routines.

::

    void shiftZ(int jx, int jy, double zangle);

``Field``
---------

The two main types are ``Field2D``, and ``Field3D``. Their main
functions are to provide an easy way to manipulate data; they take care
of all memory management, and most looping over grid-points in algebraic
expressions. The 2D field implementation is relatively simple, but more
optimisations are used in the 3D field implementation because they are
much larger (factor of :math:`\sim 100`).

To handle time-derivatives, and enable expressions to be written in the
following form:

::

    ddt(Ni) = -b0xGrad_dot_Grad(phi, Ni);

fields (and vectors, see below) have a function:

::

    Field3D* timeDeriv();

which returns a pointer to the field holding the time-derivative of this
variable. This function ensures that this field is unique using a
singleton pattern.

``Vector``
----------

Vector classes build on the field classes, just using a field to
represent each component.

To handle time-derivatives of vectors, some care is needed to ensure
that the time-derivative of each vector component points to the same
field as the corresponding component of the time-derivative of the
vector:

::

    ddt(v.x) = ddt(v).x

``dcomplex``
------------

Several parts of the BOUT++ code involve FFTs and are therefore much
easier to write using complex numbers. Unfortunately, the C++ complex
library also tries to define a ``real`` type, which is already defined
by PVODE. Several work-arounds were tried, some of which worked on some
systems, but it was easier in the end to just implement a new class
``dcomplex`` to handle complex numbers.

.. _sec-memorymanage:

Memory management
-----------------

This code has been thoroughly tested/debugged, and should only be
altered with great care, since just about every other part of BOUT++
depends on this code working correctly. Two optimisations used in the
data objects to speed up code execution are memory recycling, which
eliminates allocation and freeing of memory; and copy-on-change, which
minimises unnecessary copying of data.

Both of these optimisations are done “behind the scenes”, hidden from
the remainder of the code, and are illustrated in :numref:`fig-memory`:

.. _fig-memory:
.. figure:: ../figs/memory.*
   :alt: Memory handling in BOUT++

   Memory handling in BOUT++. Memory allocation and freeing is
   eliminated by recycling memory blocks, and assignments without
   changes (``A = B``) do not result in copying data, only pointers to
   the data. Both these optimisations are handled internally, and are
   invisible to the programmer.

The objects (A,B,C) accessed by the user in operations discussed in the
previous section act as an interface to underlying data (a,b). Memory
recycling can be used because all the scalar fields are the same size
(and vector fields are implemented as a set of 3 scalar fields). Each
class implements a global stack of available memory blocks. When an
object is assigned a value, it attempts to grab one of these memory
blocks, and if none are available then a new block is allocated. When an
object is destroyed, its memory block is not freed, but is put onto the
stack. Since the evaluation of the time-derivatives involves the same
set of operations each time, this system means that memory is only
allocated the first time the time-derivatives are calculated, after
which the same memory blocks are re-used. This eliminates the often slow
system calls needed to allocate and free memory, replacing them with
fast pointer manipulation.

Copy-on-change (reference counting) further reduces memory useage and
unnecessary copying of data. When one field is set equal to another
(e.g. ``Field3D A = B`` in :numref:`fig-memory`), no data is copied, only
the reference to the underlying data (in this case both A and B point to
data block a). Only when one of these objects is modified is a second
memory block used to store the different value. This is particularly
useful when returning objects from a routine. Usually this would involve
copying data from one object to another, and then destroying the
original copy. Using reference counting this copying is eliminated.

Global field gather / scatter
-----------------------------

In BOUT++ each processor performs calculations on a sub-set of the mesh,
and communicates with other processors primarily through exchange of
guard cells (the ``mesh->commmunicate`` function). If you need to gather
data from the entire mesh onto a single processor, then this can be done
using either 2D or 3D ``GlobalFields`` .

First include the header file

::

    #include <bout/globalfield.hxx>

which defines both ``GlobalField2D`` and ``GlobalField3D`` . To create a
3D global field, pass it the mesh pointer:

::

      GlobalField3D g3d(mesh);

By default all data will be gathered onto processor 0. To change this,
specify which processor the data should go to as the second input

::

      GlobalField3D g3d(mesh, processor);

Gather and scatter methods are defined:

::

      Field3D localData;
      // Set local data to some value

      g3d.gather(localData);  // Gathers all data onto one processor

      localData = g3d.scatter(); // Scatter data back

**Note:** Boundary guard cells are **not** handled by the scatter step,
as this would mean handling branch-cuts etc. To obtain valid data in the
guard and Y boundary cells, you will need to communicate and set Y
boundaries.

**Note:** Gather and Scatter are global operations, so all processors
must call these functions.

Once data has been gathered, it can be used on one processor. To check
if the data is available, call the method ``dataIsLocal()``, which will
return ``true`` only on one processor

::

      if(g3d.dataIsLocal()) {
        // Data is available on this processor

      }

The sizes of the global array are available through ``xSize()``,
``ySize()`` and ``zSize()`` methods. The data itself can be accessed
indirectly using ``(x,y,z)`` operators:

::

      for(int x=0; x<g3d.xSize(); x++)
        for(int y=0; y<g3d.ySize(); y++)
          for(int z=0; z<g3d.zSize(); z++)
            output.write("Value at (%d,%d,%d) is %e\n",
            x,y,z,
            g3d(x,y,z) );

or by getting a pointer to the underlying data, which is stored as a 1D
array:

::

      BoutReal *data = g3d.getData();
      nx = g3d.xSize();
      ny = g3d.ySize();
      nz = g3d.zSize();

      data[x*ny*nz + y*nz + z]; // Value at g3d(x,y,z)

See the example ``examples/test-globalfield`` for more examples.

.. _sec-iterating:

Iterating over fields
---------------------

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

.. _sec-rangeiterator:

Iterating over ranges
---------------------

The boundary of a processor’s domain may consist of a set of disjoint
ranges, so the mesh needs a clean way to tell any code which depends on
the boundary how to iterate over it. The ``RangeIterator`` class in
``include/bout/sys/range.hxx`` and ``src/sys/range.cxx`` provides this.

RangeIterator can represent a single continuous range, constructed by
passing the minimum and maximum values.

::

    RangeIterator it(1,4);  // Range includes both end points
    for(it.first(); !it.isDone(); it.next())
      cout << it.ind; // Prints 1234

A more canonical C++ style is also supported, using overloaded ``++``,
``*``, and ``!=`` operators:

::

    for(it.first(); it != RangeIterator::end(); it++)
      cout << *it; // Prints 1234

where ``it++`` is the same as ``it.next()``, and ``*it`` the same as
``it.ind``.

To iterate over several ranges, ``RangeIterator`` can be constructed
with the next range as an argument:

::

    RangeIterator it(1,4, RangeIterator(6,9));
    for(it.first(); it != RangeIterator::end(); it++)
      cout << *it; // Prints 12346789

and these can be chained together to an arbitrary depth.

To support statements like

::

    for(RangeIterator it = mesh->iterateBndryLowerY(); !it.isDone(); it++)
      ...

the initial call to ``first()`` is optional, and everything is
initialised in the constructor.

.. _sec-fieldops:

Field2D/Field3D Arithmetic Operators
------------------------------------

The arithmetic operators (``+``, ``-``, ``/``, ``*``) for ``Field2D``
and ``Field3D`` are generated automatically using the `Jinja`_
templating system. This requires Python 3 (2.7 may work, but only 3 is
supported).

Because this is fairly low-level code, and we don't expect it to
change very much, the generated code is kept in the git
repository. This has the benefit that Python and Jinja are not needed
to build BOUT++, only to change the ``Field`` operator code.

.. warning:: You should not modify the generated code
             directly. Instead, modify the template and re-generate
             the code. If you commit changes to the template and/or
             driver, make sure to re-generate the code and commit it
             as well

The Jinja template is in ``src/field/gen_fieldops.jinja``, and the
driver is ``src/field/gen_fieldops.py``. The driver loops over every
combination of ``BoutReal``, ``Field2D``, ``Field3D`` (collectively
just "fields" here) with the arithmetic operators, and uses the
template to generate the appropriate code. There is some logic in the
template to handle certain combinations of the input fields: for
example, for the binary infix operators, only check the two arguments
are on identical meshes if neither is ``BoutReal``.

To install Jinja:

.. code-block:: console

   $ pip3 install --user Jinja2

To re-generate the code, there is a ``make`` target for
``gen_fieldops.cxx`` in ``src/field/makefile``. This also tries to
apply ``clang-format`` in order to keep to a consistent code style.

.. note:: ``clang-format`` is bundled with ``clang``. This should be
          available through your system package manager. If you do not
          have sufficient privileges on your system, you can install
          it from the source `clang`_. One of the BOUT++ maintainers
          can help apply it for you too.

.. _Jinja: http://jinja.pocoo.org/
.. _clang: https://clang.llvm.org/

