#!/usr/bin/env python3

# MMS test for differential operators (that use the metric)

import boutcore

import numpy
import scipy.interpolate
from sys import exit

gridfiles = [('21712_32x32.grd.nc',32), ('21712_64x64.grd.nc',64), ('21712_128x128.grd.nc',128), ('21712_256x256.grd.nc',256)]
mxg = 2
myg = 2
plot_error = True

#boutcore.init('-q -q -q -q')
boutcore.init()

testfunc = 'sin(2*pi*x + y + z)'
boutcore_operator = boutcore.Grad_par
order = 2
stagger = None
dimensions = 'xyz'
method = None

results = []
for twistshift, paralleltransform in [('true', 'identity'), ('false', 'shifted'), ('false', 'identity'), ('true', 'shifted')]:
    print('twistshift='+twistshift+' paralleltransform='+paralleltransform)
    output_list = []
    for gridfile, nz in gridfiles:
        print(gridfile)

        # set options
        boutcore.setOption('twistshift', twistshift, force=True)
        boutcore.setOption('testmesh:paralleltransform', paralleltransform, force=True)
        boutcore.setOption('mxg', str(mxg), force=True)
        boutcore.setOption('myg', str(myg), force=True)
        boutcore.setOption('mz', str(nz), force=True)
        # set up mesh input
        boutcore.setOption('testmesh:file', 'grids/'+gridfile, force=True)
        if stagger is None:
            boutcore.setOption('testmesh:staggergrids', str('false'), force=True)
        else:
            boutcore.setOption('testmesh:staggergrids', str('true'), force=True)

        # create new Mesh object
        mesh = boutcore.Mesh(section='testmesh')

        if stagger is None:
            inloc = 'CENTRE'
            outloc = 'CENTRE'
        else:
            inloc = stagger[0]
            outloc = stagger[1]

        # calculate result of differential operator using BOUT++ implementation
        bout_input = boutcore.create3D(testfunc, mesh, outloc=inloc) - mesh.coordinates().zShift
        mesh.communicate(bout_input)

        if method is None:
            bout_result = boutcore_operator(bout_input, outloc=outloc)
        else:
            bout_result = boutcore_operator(bout_input, outloc=outloc, method=method)
        mesh.communicate(bout_result)

        output_list.append(bout_result.get()[mxg:-mxg,myg:-myg,:]) # numpy arrays of output

    lowres = output_list[0]
    nx,ny,nz = lowres.shape
    gridx = numpy.linspace(0.5,nx-0.5,nx)/nx
    gridy = numpy.linspace(0.5,ny-0.5,ny)/ny
    gridz = numpy.linspace(0.,nz-1.,nz)/nz
    for i in range(1,len(output_list)):
        output = output_list[i]
        this_nx,this_ny,this_nz = output.shape
        this_gridx = numpy.linspace(0.5,this_nx-0.5,this_nx)/this_nx
        this_gridy = numpy.linspace(0.5,this_ny-0.5,this_ny)/this_ny
        this_gridz = numpy.linspace(0.,this_nz-1.,this_nz)/this_nz

        interp_output = scipy.interpolate.Interpn([this_gridx, this_gridy, this_gridz], output, [gridx, gridy, gridz])
        output_list[i] = interp_output

    error_list = [numpy.max(numpy.abs(x-output_list[-1])) for x in output_list[:-1]]

    print(error_list)
    logerrors = numpy.log(error_list[-2]/error_list[-1])
    logspacing = numpy.log(ngrids[-1]/ngrids[-2])
    convergence = logerrors/logspacing

    if order-.1 < convergence < order+.2:
        results.append('pass')
        results.append('pass: '+str(boutcore_operator)+' is  working correctly for '+inloc+'->'+outloc+' twistshift='+twistshift+' paralleltransform='+paralleltransform+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.')
        #results.append('Proc #'+str(mesh.getYProcIndex())+' --- pass: '+str(boutcore_operator)+' is  working correctly for '+inloc+'->'+outloc+' twistshift='+twistshift+' paralleltransform='+paralleltransform+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.')
    else:
        if plot_error:
            from matplotlib import pyplot
            pyplot.loglog(1./ngrids, error_list)
            pyplot.show()
            from boututils.showdata import showdata
            showdata(error)
        results.append(str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' twistshift='+twistshift+' paralleltransform='+paralleltransform+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.')
        #results.append('Proc #'+str(mesh.getYProcIndex())+' --- '+str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' twistshift='+twistshift+' paralleltransform='+paralleltransform+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.')

# check results of tests
fail = False
for result in results:
    if result is not 'pass':
        print(result)
        fail = True
mesh.communicate(bout_result) # use this like MPI_Barrier() to synchronise processes
if fail:
    exit(1)
else:
    print('pass')
    exit(0)
