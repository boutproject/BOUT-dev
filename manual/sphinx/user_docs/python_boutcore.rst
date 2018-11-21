The python boutcore module
==========================

Installing
----------

Installing boutcore can be tricky.
Ideally it should be just

.. code-block:: bash

   ./configure --enable-shared
   make -j 4 python


but getting all the
dependencies can be difficult.
``make python`` creates the python3 module.

If problems arise, it might be worth checking a copy of the bout
module out, to reduce the risk of causing issues with the old bout
installation. This is especially true if you are trying to run
boutcore not on compute nodes of a super computer but rather on
post-processing/login/... nodes.

To use boutcore on the login node, a self compiled version of mpi may be
required, as the provided one may be only for the compute nodes.
Further, numpy header files are required, therefore numpy needs to be
compiled as well.
Further, the header files need to be exposed to the boutcore cython
compilation, e.g. by adding them to ``_boutcore_build/boutcore.pyx.in``.
It seems both ``NUMPY/numpy/core/include`` and
``NUMPY/build/src.linux-x86_64-2.7/numpy/core/include/numpy`` need to be
added, where ``NUMPY`` is the path of the numpy directory.
For running boutcore on the post processing nodes, fftw3 needs to be
compiled as well, if certain fftw routines are used. Note, fftw needs
to be configured with ``--enable-shared``.

After installing mpi e.g. in ``~/local/mpich``, bout needs to be
configured with something like:
``./configure --enable-shared MPICC=~/local/mpich/bin/mpicc MPICXX=~/local/mpich/bin/mpicxx --with-fftw=~/local/fftw/``

``--enable-shared`` is required, so that pvode etc. is compiles as position
independent code.

If you are running fedora - you can install pre-build binaries:

.. code-block:: bash

   sudo dnf copr enable davidsch/bout
   sudo dnf install python3-bout++-mpich
   module load mpi/mpich-$(arch)


Purpose
-------

The boutcore module exposes (part) of the BOUT++ C++ library to python.
It allows to calculate e.g. BOUT++ derivatives in python.

State
-----
Field3D and Field2D are working. If other fields are needed, please open an issue.
Fields can be accessed directly using the [] operators, and give a list of slice objects.
The get all data, ``f3d.getAll()`` is equivalent to ``f3d[:,:,]`` and returns a numpy array.
This array can be addressed with
e.g. ``[]`` operators, and then the field can be set again with
``f3d.setAll(numpyarray)``.
It is also possible to set a part of an Field3D with the ``[]`` operators.
Addition, multiplication etc. are all available.
The derivatives should all be working, if find a missing one, please open an issue.
Vectors are not exposed yet.

Functions
---------

.. automodule:: boutcore
   :members:


Examples
--------
Some trivial post processing:

.. code-block:: python

   import boutcore
   import numpy as np
   args="-d data -f BOUT.settings -o BOUT.post".split(" ")
   boutcore.init(args)
   dens=boutcore.Field3D.fromCollect("n",path="data")
   temp=boutcore.Field3D.fromCollect("T",path="data")
   pres=dens*temp
   dpdz=boutcore.DDZ(pres,outloc="CELL_ZLOW")



A simple MMS test:

.. code-block:: python

   import boutcore
   import numpy as np
   boutcore.init("-d data -f BOUT.settings -o BOUT.post")
   for nz in [64,128,256]:
       boutcore.setOption("meshz:nz","%d"%nz)
       mesh=boutcore.Mesh(OptionSection="meshz")
       f=boutcore.create3D("sin(z)",mesh)
       sim=boutcore.DDZ(f)
       ana=boutcore.create3D("cos(z)",mesh)
       err=sim-ana
       err=boutcore.max(boutcore.abs(err))
       errors.append(err)


A real example - unstagger data:

.. code-block:: python

   import boutcore
   boutcore.init("-d data -f BOUT.settings -o BOUT.post")
   # uses location from dump - is already staggered
   upar=boutcore.Field3D.fromCollect("Upar")
   upar=boutcore.interp_to(upar,"CELL_CENTRE")
   # convert to numpy array
   upar=upar.getAll()


A real example - check derivative contributions:

.. code-block:: python

   #!/usr/bin/env python

   from boutcore import *
   import numpy as np
   from netCDF4 import Dataset
   import sys

   if len(sys.argv)> 1:
       path=sys.argv[1]
   else:
       path="data"

   times=collect("t_array",path=path)

   boutcore.init("-d data -f BOUT.settings -o BOUT.post")
   with Dataset(path+'/vort.nc', 'w', format='NETCDF4') as outdmp:
      phiSolver=Laplacian()
      phi=Field3D.fromCollect("n",path=path,tind=0,info=False)
      zeros=phi.getAll()*0
      phi.setAll(zeros)
      outdmp.createDimension('x',zeros.shape[0])
      outdmp.createDimension('y',zeros.shape[1])
      outdmp.createDimension('z',zeros.shape[2])
      outdmp.createDimension('t',None)
      t_array_=outdmp.createVariable('t_array','f4',('t'))
      t_array_[:]=times
      ExB     = outdmp.createVariable('ExB'    ,'f4',('t','x','y','z'))
      par_adv = outdmp.createVariable('par_adv','f4',('t','x','y','z'))
      def setXGuards(phi,phi_arr):
          for z in range(tmp.shape[2]):
              phi[0,:,z]=phi_arr
              phi[1,:,z]=phi_arr
              phi[-2,:,z]=phi_arr
              phi[-1,:,z]=phi_arr
      with open(path+"/equilibrium/phi_eq.dat","rb") as inf:
          phi_arr=np.fromfile(inf,dtype=np.double)
          bm="BRACKET_ARAKAWA_OLD"

          for tind in range(len(times)):
              vort     = Field3D.fromCollect("vort"     ,path=path,tind=tind,info=False)
              U        = Field3D.fromCollect("U"        ,path=path,tind=tind,info=False)
              setXGuards(phi,phi_arr)
              phi=phiSolver.solve(vort,phi)
              ExB[tind,:,:,:]=(-bracket(phi, vort, bm, "CELL_CENTRE")).getAll()
              par_adv[tind,:,:,:]=(- Vpar_Grad_par(U, vort)).getAll()



Functions - undocumented
------------------------

.. automodule:: boutcore
   :undoc-members:

Functions - special and inherited
---------------------------------

.. automodule:: boutcore
   :special-members:
   :inherited-members:
