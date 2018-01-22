Miscellaneous
=============

Other small modules which don’t really fit into any system, but are
needed.

Printing messages
-----------------

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
