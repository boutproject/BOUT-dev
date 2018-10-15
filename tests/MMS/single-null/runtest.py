#!/usr/bin/env python3

# MMS test for differential operators (that use the metric)

import boutcore

from boutdata.mms_alternate import *

import numpy
import sympy
from sys import exit

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--short', action='store_true', default=False)
args = parser.parse_args()
full_test = not args.short

min_exponent = 6
max_exponent = 7
ngrids = numpy.logspace(min_exponent, max_exponent, num=max_exponent-min_exponent+1, base=2, dtype=int)
default_n = 4
mxg = 2
myg = 2
plot_error = False
test_derivs = full_test

boutcore.init('-q -q -q -q')

def test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore_operator, symbolic_operator, order, method=None):
    """
    Test a set of parameters

    Parameters
    ----------

    twistshift : str
        'true' or 'false' value for twistshift option
    paralleltransform : str
        value for paralleltransform option
    dimensions : str
        string containing any of 'x', 'y' and 'z': the grid must be refined in
        the directions contained in this argument in order to converge. E.g.
        Grad_par requires refinement only in the y-direction.
    stagger : tuple of str, None
        inloc/outloc for the operator, or None for staggergrids=false
    ngrids : numpy.array(int)
        Array of grid sizes to use. All directions being refined are given the
        same size.
    testfunc : sympy expression
        The input that will be given to the function being tested
    boutcore_operator : function
        function from boutcore to be tested
    symbolic_operator : function
        function using sympy to do the symbolic equivalent of boutcore_operator
    order : int
        expected order of convergence of boutcore_operator
    method : str or None
        method argument for operator, e.g. 'C2'. None to use value from
        BOUT.inp (or default if no option set).
    """

    print('twistshift='+twistshift+' paralleltransform='+paralleltransform)
    print(boutcore_operator)

    # geometry for simple circular tokamak
    if paralleltransform == 'identity':
        tokamak = SimpleTokamak(psiN0=.95, lower_legs_fraction=0.25)
    elif paralleltransform == 'shifted':
        tokamak = SimpleTokamak(psiN0=.95, lower_legs_fraction=0.25, shifted=True)
    else:
        raise ValueError('Unhandled paralleltransform='+paralleltransform)
    # rescale x and y coordinates so dx and dy are not constants
    tokamak.set_scalex(1 + .1*sin(2*pi*metric.x+metric.y))
    tokamak.set_scaley(1 + .1*sin(2*pi*metric.x-metric.y))
    # re-calculate metric terms
    tokamak.metric()

    error_list = []
    for n in ngrids:
        print('n =',n)

        if 'x' in dimensions:
            nx = n
        else:
            nx = default_n
        if 'y' in dimensions:
            ny = n
        else:
            ny = default_n
        if 'z' in dimensions:
            nz = n
        else:
            nz = default_n

        tokamak.setNs(nx=nx, ny=ny, nz=nz, mxg=mxg)

        # set options
        boutcore.setOption('twistshift', twistshift, force=True)
        boutcore.setOption('testmesh:paralleltransform', paralleltransform, force=True)
        boutcore.setOption('mxg', str(mxg), force=True)
        boutcore.setOption('myg', str(myg), force=True)
        # set up mesh input
        boutcore.setOption('testmesh:nx', exprToStr(nx+2*mxg), force=True)
        boutcore.setOption('testmesh:ny', exprToStr(ny), force=True)
        boutcore.setOption('testmesh:nz', exprToStr(nz), force=True)
        boutcore.setOption('testmesh:dx', exprToStr(tokamak.dx), force=True)
        boutcore.setOption('testmesh:dy', exprToStr(tokamak.dy), force=True)
        boutcore.setOption('testmesh:dz', exprToStr(tokamak.dz), force=True)
        boutcore.setOption('testmesh:g11', exprToStr(metric.g11), force=True)
        boutcore.setOption('testmesh:g22', exprToStr(metric.g22), force=True)
        boutcore.setOption('testmesh:g33', exprToStr(metric.g33), force=True)
        boutcore.setOption('testmesh:g12', exprToStr(metric.g12), force=True)
        boutcore.setOption('testmesh:g13', exprToStr(metric.g13), force=True)
        boutcore.setOption('testmesh:g23', exprToStr(metric.g23), force=True)
        boutcore.setOption('testmesh:g_11', exprToStr(metric.g_11), force=True)
        boutcore.setOption('testmesh:g_22', exprToStr(metric.g_22), force=True)
        boutcore.setOption('testmesh:g_33', exprToStr(metric.g_33), force=True)
        boutcore.setOption('testmesh:g_12', exprToStr(metric.g_12), force=True)
        boutcore.setOption('testmesh:g_13', exprToStr(metric.g_13), force=True)
        boutcore.setOption('testmesh:g_23', exprToStr(metric.g_23), force=True)
        boutcore.setOption('testmesh:zShift', exprToStr(metric.zShift), force=True)
        boutcore.setOption('testmesh:shiftAngle', exprToStr(metric.shiftAngle), force=True)
        boutcore.setOption('testmesh:ixseps1', exprToStr(tokamak.ixseps1), force=True)
        boutcore.setOption('testmesh:ixseps2', exprToStr(tokamak.ixseps2), force=True)
        boutcore.setOption('testmesh:jyseps1_1', exprToStr(tokamak.jyseps1_1), force=True)
        boutcore.setOption('testmesh:jyseps1_2', exprToStr(tokamak.jyseps1_2), force=True)
        boutcore.setOption('testmesh:jyseps2_1', exprToStr(tokamak.jyseps2_1), force=True)
        boutcore.setOption('testmesh:jyseps2_2', exprToStr(tokamak.jyseps2_2), force=True)
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
        bout_input = boutcore.create3D(exprToStr(testfunc.subs(metric.z, metric.z + metric.zShift)), mesh, outloc=inloc)
        mesh.communicate(bout_input)

        if method is None:
            bout_result = boutcore_operator(bout_input, outloc=outloc)
        else:
            bout_result = boutcore_operator(bout_input, outloc=outloc, method=method)
        mesh.communicate(bout_result)
        # calculate result of differential operator symbolically, then convert to boutcore.Field3D/Field2D
        if paralleltransform == 'identity':
            analytic_func = symbolic_operator(testfunc.subs(metric.z, metric.z + metric.zShift))
        elif paralleltransform == 'shifted':
            analytic_func = symbolic_operator(testfunc).subs(metric.z, metric.z + metric.zShift)
        analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh, outloc=outloc)

        # calculate max error
        error = bout_result - analytic_result # as Field3D/Field2D
        error = error.get()[mxg:-mxg, myg:-myg] # numpy array, without guard cells
        error_list.append(numpy.max(numpy.abs(error))) # max error

    print(error_list)
    logerrors = numpy.log(error_list[-2]/error_list[-1])
    logspacing = numpy.log(ngrids[-1]/ngrids[-2])
    convergence = logerrors/logspacing

    if order-.1 < convergence < order+.2:
        #return['pass']
        #return['pass: '+str(boutcore_operator)+' is  working correctly for '+inloc+'->'+outloc+' twistshift='+twistshift+' paralleltransform='+paralleltransform+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.']
        return ['Proc #'+str(mesh.getYProcIndex())+' --- pass: '+str(boutcore_operator)+' is  working correctly for '+inloc+'->'+outloc+' twistshift='+twistshift+' paralleltransform='+paralleltransform+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.']
    else:
        if plot_error:
            yproc = mesh.getYProcIndex()
            from matplotlib import pyplot
            pyplot.loglog(1./ngrids, error_list)
            pyplot.title("proc = "+str(yproc))
            pyplot.show()
            from boututils.showdata import showdata
            showdata(error, titles=["proc = "+str(yproc)])
            #showdata([bout_result.get()[mxg:-mxg, myg:-myg], analytic_result.get()[mxg:-mxg, myg:-myg], error], titles=["proc = "+str(yproc)+" bout", "proc = "+str(yproc)+" analytic", "proc = "+str(yproc)+" error"])
        #return[str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' twistshift='+twistshift+' paralleltransform='+paralleltransform+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.']
        return['Proc #'+str(mesh.getYProcIndex())+' --- '+str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' twistshift='+twistshift+' paralleltransform='+paralleltransform+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.']

#testfunc = sin(2*pi*metric.x + metric.y + metric.z)
testfunc = sin(z)
order = 2
dimensions = 'xyz'
stagger = None

results = []

#for twistshift, paralleltransform in [('true', 'identity'), ('false', 'shifted')]:
for twistshift, paralleltransform in [('false', 'shifted')]:
    results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.Grad_par, Grad_par, order)
    if test_derivs:
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.DDX, DDX, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.DDY, DDY, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.DDZ, DDZ, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DX2, D2DX2, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DY2, D2DY2, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DZ2, D2DZ2, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DXDY, D2DXDY, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DYDX, D2DYDX, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DXDZ, D2DXDZ, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DZDX, D2DZDX, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DYDZ, D2DYDZ, order)
        results += test_operator(twistshift, paralleltransform, dimensions, stagger, ngrids, testfunc, boutcore.D2DZDY, D2DZDY, order)

# check results of tests
fail = False
for result in results:
    if result is not 'pass':
        print(result)
        fail = True
boutcore.MPI_Barrier()
if fail:
    exit(1)
else:
    print('pass')
    exit(0)
