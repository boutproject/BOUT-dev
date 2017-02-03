.. _sec-3to4:

Updating Physics Models from v3 to v4
=====================================

Version 4.0.0 of BOUT++ introduced several features which break backwards
compatibility. If you already have physics models, you will most likely need to
update them to work with version 4. The main breaking changes which you are
likely to come across are:

* Using round brackets ``()`` instead of square brackets ``[]`` for indexing
  fields

* Moving components of :cpp:class:`Mesh` related to the metric tensor and "real
  space" out into a new object, :cpp:class:`Coordinates`

* Changed some :cpp:class:`Field3D` member functions into non-member functions

A new tool is provided, ``bin/bout_3to4.py``, which can identify these changes,
and fix most of them automatically. Simply run this program on your physic model
to see how to update it to work with version 4:

.. code-block:: bash

   $ ${BOUT_TOP}/bin/bout_3to4.py my_model.cxx

The output of this command will show you how to fix each problem it
identifies. To automatically apply them, you can use the ``--replace`` option:

.. code-block:: bash

   $ ${BOUT_TOP}/bin/bout_3to4.py --replace my_model.cxx

Also in version 4 is a new syntax for looping over each point in a field. See
:ref:`sec-iterating` for more information.
