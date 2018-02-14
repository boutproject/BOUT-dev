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
          it from the source `here`_. One of the BOUT++ maintainers
          can help apply it for you too.

.. _Jinja: http://jinja.pocoo.org/
.. _here: https://clang.llvm.org/
