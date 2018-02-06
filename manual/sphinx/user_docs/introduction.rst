.. _sec-userguide:

Introduction
============

BOUT++ is a C++ framework for writing plasma fluid simulations with an
arbitrary number of equations in 3D curvilinear coordinates . It has
been developed from the original **BOU**\ ndary **T**\ urbulence 3D
2-fluid edge simulation code written by X.Xu and M.Umansky at LLNL.

Though designed to simulate tokamak edge plasmas, the methods used are
very general and almost any metric tensor can be specified, allowing the
code to be used to simulate (for example) plasmas in slab, sheared slab,
and cylindrical coordinates. The restrictions on the simulation domain
are that the equilibrium must be axisymmetric (in the z coordinate), and
that the parallelisation is done in the :math:`x` and :math:`y`
(parallel to :math:`\mathbf{B}`) directions.

The aim of BOUT++ is to automate the common tasks needed for simulation
codes, and to separate the complicated (and error-prone) details such as
differential geometry, parallel communication, and file input/output
from the user-specified equations to be solved. Thus the equations being
solved are made clear, and can be easily changed with only minimal
knowledge of the inner workings of the code. As far as possible, this
allows the user to concentrate on the physics, rather than worrying
about the numerics. This doesn’t mean that users don’t have to think
about numerical methods, and so selecting differencing schemes and
boundary conditions is discussed in this manual. The generality of the
BOUT++ of course also comes with a limitation: although there is a large
class of problems which can be tackled by this code, there are many more
problems which require a more specialised solver and which BOUT++ will
not be able to handle. Hopefully this manual will enable you to test
whether BOUT++ is suitable for your problem as quickly and painlessly as
possible.

This manual is written for the user who wants to run (or modify)
existing plasma models, or specify a new problem (grid and equations) to
be solved. In either case, it’s assumed that the user isn’t all that
interested in the details of the code. For a more detailed descriptions
of the code internals, see the developer and reference guides. After
describing how to install BOUT++ (section :ref:`sec-install`), run the test
suite (section :ref:`sec-runtestsuite`) and a few examples
(section :ref:`sec-running`, more detail in section :ref:`sec-examples`),
increasingly sophisticated ways to modify the problem being solved are
introduced. The simplest way to modify a simulation case is by altering
the input options, described in section :ref:`sec-options`. Checking that the
options are doing what you think they should be by looking at the output
logs is described in section :ref:`sec-running`, and an overview of the IDL
analysis routines for data post-processing and visualisation is given in
section :ref:`sec-output`. Generating new grid files, particularly for
tokamak equilibria, is described in section :ref:`sec-gridgen`.

Up to this point, little programming experience has been assumed, but
performing more drastic alterations to the physics model requires
modifying C++ code. Section :ref:`sec-equations` describes how to write a new
physics model specifying the equations to be solved, using ideal MHD as
an example. The remaining sections describe in more detail aspects of
using BOUT++: section :ref:`sec-diffops` describes the differential operators
and methods available; section :ref:`sec-staggergrids` covers the
experimental staggered grid system.

Various sources of documentation are:

-  Most directories in the BOUT++ distribution contain a README file.
   This should describe briefly what the contents of the directory are
   and how to use them.

-  This user’s manual, which goes through BOUT++ from a user’s point of
   view

-  The developer’s manual, which gives details of the internal working
   of the code.

-  The reference guide, which summarises functions, settings etc.
   Intended more for quick reference rather than a guide.

- Most of the code contains Doxygen comment tags (which are slowly
  getting better). Running `doxygen <www.doxygen.org>`_ on these files
  should therefore generate an HTML reference. This is probably going
  to be the most up-to-date documentation.

License and terms of use
------------------------

Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu

BOUT++ is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BOUT++ is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with BOUT++.  If not, see
<https://www.gnu.org/licenses/>.

A copy of the LGPL license is in COPYING.LESSER. Since this is based
on (and refers to) the GPL, this is included in COPYING.

BOUT++ is free software, but since it is a scientific code we also ask
that you show professional courtesy when using this code:

#. Since you are benefiting from work on BOUT++, we ask that you submit
   any improvements you make to the code to us by emailing Ben Dudson at
   bd512@york.ac.uk

#. If you use BOUT++ results in a paper or professional publication, we
   ask that you send your results to one of the BOUT++ authors first so
   that we can check them. It is understood that in most cases if one or
   more of the BOUT++ team are involved in preparing results then they
   should appear as co-authors.

#. Publications or figures made with the BOUT++ code should
   acknowledge the BOUT++ code by citing `B.Dudson
   et. al. Comp.Phys.Comm 2009`_ and/or other BOUT++ papers. See the
   file CITATION for details.

..
   .. toctree::
      :maxdepth: 1

      overview
      getting_started
      advanced_install
      running_bout
      makefiles
      output_and_post
      bout_options
      variable_init
      time_integration
      boundary_options
      iterating
      input_grids
      laplacian
      fluid_equations
      fluid_equations_2
      object_orientated_interface
      differential_operators
      staggered_grids
      advanced_methods
      eigenvalue_solver
      testing
      notes
      machine_install
      aix
      bout_functions_for_physics
      idl
      python
      fourier_transform_derivatives
      examples


.. _B.Dudson et. al. Comp.Phys.Comm 2009: https://www.sciencedirect.com/science/article/B6TJ5-4VTCM95-3/2/ed200cd23916d02f86fda4ce6887d798
