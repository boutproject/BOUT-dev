Contributing to BOUT++
======================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

If you would like to help contribute to BOUT++ then there are many
things you can do which will make a difference. There are projects large
and small which a student or researcher could use to get started with
BOUT++ and get more familiar with the code. You don’t need to be
particularly familar with BOUT++ or C++ to work on many of these. You
can see a current list of outstanding bugs and feature requests on the
`GitHub issue page <https://github.com/boutproject/BOUT-dev/issues>`__.

This and following sections describe the core of BOUT++, and are intended
for anyone who wants to work on improving BOUT++. For a
general introduction, and instructions for using BOUT++ see the :ref:`User guide
<sec-userguide>`. The user’s guide assumes only minimal knowledge of C++, and
provides only those details needed to use BOUT++.

We use doxygen comments to document the code in the source files. This is then
built using breathe_ and sphinx_ into the document you're reading now. The
API documentation is in :ref:`api-ref`.


.. _breathe: https://breathe.readthedocs.io
.. _sphinx: http://www.sphinx-doc.org/

House rules
-----------

As production codes go, BOUT++ is not particularly big,
but it is definitely large enough that keeping the code ‘clean’ and
understandable is necessary. This is vital if many people are going to
work on the code, and also greatly helps code debugging and
verification. There are therefore a few house rules to keep in mind when
modifying the BOUT++ code.

When modifying the core BOUT++ code, please keep in mind that this
portion of the code is intended to be general (i.e. independent of any
particular physical system of equations), and to be used by a wide range
of users. Making code clear is also more important in this section than
the physics model since the number of developers is potentially much
greater.

Here are some rules for editing the core BOUT++ code:

-  **NO FORTRAN**. EVER. Though it may be tempting for scientific
   programmers to use a little Fortran now and then, please please don’t
   put any into BOUT++. Use of Fortran, particularly when mixed with
   C/C++, is the cause of many problems in porting and modifying codes.

-  If a feature is needed to study a particular system, only include it
   in the core code if it is more generally applicable, or cannot be put
   into the physics module.

-  If you add a new feature, function, class member, etc. you must also
   include doxygen comments that explain what each new thing does.
   Similarly, if a change you make would affect e.g. a function’s
   arguments, please ensure that you keep the documentation up-to-date
   with the code. See the section on coding style for best practices in
   this regard. If you submit a pull request that doesn’t add or update
   documentation where appropriate, we may ask you to do so before it is
   merged.

-  As well as documentation for new features, you must also include a
   representative test to ensure that it works correctly. Please see the
   `tests README <./tests/README.md>`__ for more information on tests in
   BOUT++. Prefer to write unit tests that check the feature at the
   function level, rather than integrated tests that require setting up
   a whole physics model.

Development workflow using Git
------------------------------

The workflow we use is essentially
`“gitflow” <https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow>`__.

-  **master** should always be stable
-  **next** contains bleeding-edge features

All work should be done in feature branches, branched off *next*. When
complete, a pull request can be submitted.

At irregular intervals, we will create a **release** branch. No new
features go into the *release* branch - only bug fixes and
documentation.

1. Create a new branch
2. (Optional) Push it to Github to share and for backup
3. Make changes, commits
4. Submit a pull request into **next** using Github’s `Pull
   Requests <https://github.com/boutproject/BOUT-dev/pulls>`__ system

Creating a feature branch
~~~~~~~~~~~~~~~~~~~~~~~~~

First get a copy of the
`BOUT-dev <https://github.com/boutproject/BOUT-dev>`__ repository (or
git pull to update an existing copy):

.. code-block:: console

    git clone git@github.com:boutproject/BOUT-dev.git
    cd BOUT-dev

Create a new branch **myfeature**, branching from **next**. Choose a
descriptive name for **myfeature**, anything except “master” or “next”.

.. code-block:: console

    git checkout next
    git pull
    git checkout -b myfeature     # Switched to a new branch "myfeature"

Pushing to Github
~~~~~~~~~~~~~~~~~

Create a fork on Github following the instructions
`here <https://help.github.com/articles/fork-a-repo/>`__.

If you want to push your branch to BOUT-dev to share with other
developers, run::

    git push -u yourfork myfeature

This command pushes **myfeature** to your fork (named **yourfork**) of
the BOUT-dev repository, and the -u flag adds it as a remote tracking
branch. After setting up the tracking branch, you can call “git push”
without any parameters to push updates to **myfeature**.

If another developer wants to try out this branch, they will first need
to add your repository as a new remote::

    git remote add yourfork https://github.com/YourUsername/BOUT-dev.git

then they will be able to checkout your branch::

    git checkout -b myfeature yourfork/myfeature

*Note*: If you have write access to the central BOUT-dev repository, you
can push your branches there.

Making changes, commits
~~~~~~~~~~~~~~~~~~~~~~~

Now you would make changes, commit changes and push as usual:

.. code-block:: console

    ... make changes ...
    git add <files>
    git commit
    git push   # Pushes to origin/myfeature

You can switch between branches using *checkout*:

.. code-block:: console

    git checkout master    # Switch to "master"
    git checkout myfeature # Switch to "myfeature"

Merging into **next**
~~~~~~~~~~~~~~~~~~~~~

Once your feature is complete, ask other developers to have a look by
creating a `Pull
Request <https://github.com/boutproject/BOUT-dev/pulls>`__ on the
`BOUT-dev <https://github.com/boutproject/BOUT-dev>`__ page. One of the
maintainers will review your code and merge it into **next**. They may
give you comments to improve the code. You can make additional changes
and push them to the same feature branch and they will be automatically
added to the pull request.

Running Tests
~~~~~~~~~~~~~

We run many tests and checks automatically on GitHub (collectively called "CI")
for a variety of build configurations. See :ref:`sec-runtestsuite` for how to
run them locally. Running the full test suite can take some time, but it's a
good idea to run at least the unit tests locally when you're making changes:

.. code-block:: console

    cmake --build build --target check-unit-tests

Along with the automated tests, we also run things like formatters and static
analysis in CI, which may provide further feedback (see :ref:`sec-coding-style`
for more details). If your PR fails the format check on CI, please install the
formatters and run them according to the next section. Usually this is just:

.. code-block:: console

    git clang-format next

which will format just the parts of C++ files that have changed since ``next``.

Formatting and Linters
~~~~~~~~~~~~~~~~~~~~~~

For frequent developers, we strongly recommend installing the formatting and
linting tools locally and building them into your regular workflow. You can
install the majority of our developer tools at once using `uv
<https://docs.astral.sh/uv/>`_ from the top of your BOUT++ directory:

.. code-block:: console

    uv sync --only-dev --inexact

This will install all the developer tools into a virtual environment (creating
one if necessary, typically under ``.venv``). You can then activate the virtual
environment and use the tools manually, or via integrations in your editor.

We *strongly* recommend using ``uv`` to install the developer tools, as this
will ensure that you get the *exact* versions used in CI and by other
developers, but you could also use ``pip install --group dev`` if you wish.

The quickest way to run all the formatters on a PR at once is typically:

.. code-block:: console

    uv tool run prek run --from-ref next

This will run all our formatters on files that have changed since ``next``.

The tools can be updated with:

.. code-block:: console

    uv lock -U
    uv tool run sync-with-uv


Install Pre-Commit Hooks (Optional but recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We also have a `prek <https://prek.j178.dev>`_ config that can run some of
these automatically using pre-commit hooks -- that is, when you run ``git
commit``, these tools will run and abort the commit if they change any
files. You then just have to stage the new changes and try again.

Install our developer tools with ``uv`` as above, and then install the
pre-commit hook:

.. code-block:: console

    prek install

That's it! ``prek`` will install the tools into their own separate virtual
environment so they won't interfere with any you may have, and ``git`` will take
care of running them automatically.

.. _sec-coding-style:

Coding Style
------------

Code is read an order of magnitude more times than it is written. It’s
also written for *people* and not for the computer! For these reasons,
it’s important that we stick to some form of coding standards. The
following coding style guidelines broadly follow the `LLVM Coding
Standards <http://llvm.org/docs/CodingStandards.html>`__. The LLVM
Coding Standards go into more depth, and explain the reasoning behind
the guidelines more thoroughly than here. If you just follow the
guidelines below, you won’t go far wrong though.

These guidelines are intended to make code easier to read and therefore
easier to understand. Being consistent in coding style also helps
comprehension by reducing cognitive load.

We use various tools to enforce code style and quality:

- `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ for formatting
  of C++ code
- `clang-tidy <https://clang.llvm.org/extra/clang-tidy/index.html>`_ for linting
  and static analysis of C++ code
- `ruff <https://astral.sh/ruff>`_ for formatting and linting of Python code
- `cmake-format <https://github.com/cheshirekow/cmake_format>`_ for formatting
  CMake files
- `sphinx-lint <https://github.com/sphinx-contrib/sphinx-lint>`_ for linting
  documentation pages

You can install all of these tools using:

.. code-block:: console

    pip install -r requirements_dev.txt


Comments
~~~~~~~~

Comments in the code are vital to helping understanding. Comments that
are embedded in the code should explain **why** something is done,
rather than **how**.

For documenting what functions and classes do, we use
`Doxygen <http://www.doxygen.org>`__.

-  Prefer C++ style comments ``//`` over C style ``/* */``

-  Doxygen comments: use ``///``.

   Doxygen done right::

       /// Foo the bar
       ///
       /// Apply the standard foo method to @p bar
       ///
       /// Typical usage:
       ///
       ///    foo(bar, "simple", result);
       ///
       ///
       /// @param[in] quux Number of times to foo
       /// @param[out] result filled with quux fooed bars
       ///
       /// @returns true on success
       bool applyFoo(BoutReal bar, int quux, std::vector<int> &result);

-  The header files are essentially the “public” API, so prefer to put
   Doxygen comments there, rather than in the implementation. “Private”
   functions, etc., can be documented in the implementation.

Naming
~~~~~~

Naming things correctly is super important! It is also one of the
trickiest parts of coding.

Names should be *descriptive*. Code is read an order of magnitude more
often than it is written, so it is vital that it is easy to comprehend.

There are some conventions you should follow when naming things:

-  Type or class names should be nouns and be PascalCase –
   e.g. BoutReal, BoutMesh, Laplacian
-  Variable names should be nouns and snake_case – e.g. generator,
   forward_map, extra_yguards_lower
-  Functions should be verbs (i.e. actions) and camelCase – e.g. solve,
   getSection, parseString

Prefer a longer descriptive name over a shorter abbreviated one:
``inner_boundary_flags`` rather than ``inbndflgs``, ``generator`` rather
than ``gen``. It’s much easier to read and comprehend than the
abbreviated form.

Details
~~~~~~~

-  Use spaces instead of tabs. Tabs may be interpreted differently by
   different editors, making the code look badly indented and difficult
   to read. The easiest solution is just use spaces everything instead.
-  Two spaces for indentation
-  Spaces after ``if``, ``for``, etc.

   Wrong::

       if(expr){
         doSomething();
       }else{
         doOtherThing();
       }

   Right::

       if (expr) {
         doSomething();
       } else {
         doOtherThing();
       }

   This especially helps readability, making conditional statements
   stand out over function calls.
-  Braces on same line as statement:

   Wrong::

       void doFoo(bool expr)
       {
         if (expr)
         {
           doSomething();
         }
         else
         {
           doOtherThing();
         }
       }

   Right::

       void doFoo(bool expr) {
         if (expr) {
           doSomething();
         } else {
           doOtherThing();
         }
       }

   This one is more style than readability - it’s the style that the
   majority of BOUT++ already uses.
