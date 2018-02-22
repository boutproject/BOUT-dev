Data types
==========

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

