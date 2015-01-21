BOUT++ 
=======

A modular fluid simulation code written in C++. 
Aims to be able to solve a wide variety of fluid models in 
almost any curvilinear coordinate system.

In this directory
-----------------

* ABOUT                 More details and terms of use.
* CITATION              Contains the paper citation for BOUT++
* COPYING               GPL license
* COPYING.LESSER        LGPL license
* src                   The main code directory
* manual                Contains documentation, in particular the user manual.
                        Read this if you're new to BOUT / BOUT++.
* examples              Example models and test codes
* tools                 Tools for helping with analysis, mesh generation, and data managment

  * archiving           Routines for managing input/output files
                        e.g. compressing data, converting formats, and managing runs
  * idllib              Analysis codes in IDL. Add this to your IDL_PATH
                        environment variable.
  * pdb2idl             Library to read Portable Data Binary (PDB) files into IDL
  * pylib               Analysis codes in Python
  
    * boutdata        Routines to simplify accessing BOUT++ output
    * boututils       Some useful routines for accessing and plotting data

  * tokamak_grids     Code to generate input grids for tokamak equilibria
  
    * gridgen         Grid generator in IDL. Hypnotoad GUI for
                      converting G-EQDSK files into a flux-aligned
                      orthogonal grid.
    * elite           Convert ELITE .eqin files into an intermediate binary file
    * gato            Convert DSKGATO files into intermediate binary format 
    * all             Convert the intermediate binary file into BOUT++ input grid
    * shifted_circle  Produce shifted cirle equilibria input grids


License
-------

BOUT++ is released under the Lesser General Public License (LGPL). See the LICENSE file for details.

Terms of use
------------

BOUT++ is released under the LGPL, but since BOUT++ is a
scientific code we also ask that you show professional courtesy
when using this code:

1. Since you are benefiting from work on BOUT++, we ask that you
   submit any improvements you make to the code to us by emailing 
   Ben Dudson at <bd512@york.ac.uk>
2. If you use BOUT++ results in a paper or professional publication,
   we ask that you send your results to one of the BOUT++ authors
   first so that we can check them. It is understood that in most cases
   if one or more of the BOUT++ team are involved in preparing results
   then they should appear as co-authors.
3. Publications or figures made with the BOUT++ code should acknowledge the
   BOUT++ code by citing B.Dudson et. al. Comp.Phys.Comm 2009 and/or
   other BOUT++ papers. See the file CITATION for details.
