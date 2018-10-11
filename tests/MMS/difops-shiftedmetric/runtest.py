#!/usr/bin/env python3

# MMS test for differential operators (that use the metric), using the
# ShiftedMetric ParallelTransform.
# ShiftedMetric is not compatibile with y-staggering, so just test
# staggergrids=false case.
# Also operators with no y-component don't use ShiftedMetric functionality so
# skip them here.

import boutcore
from boutdata.mms_alternate import *

import numpy
import sympy
from copy import copy
from sys import exit

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--short', action='store_true', default=False)
args = parser.parse_args()
full_test = not args.short

# geometry for simple circular tokamak
tokamak = SimpleTokamak(shifted=True)

# rescale x and y coordinates so dx and dy are not constants
tokamak.set_scalex(1 + .1*sin(2*pi*metric.x+metric.y))
tokamak.set_scaley(1 + .1*sin(2*pi*metric.x-metric.y))
# re-calculate metric terms
tokamak.metric()

def test_operator(dimensions, ngrids, testfunc, boutcore_operator, symbolic_operator, order, method=None):
    """
    Test a set of parameters

    Parameters
    ----------

    dimensions : str
        String containing any of 'x', 'y' and 'z': the grid must be refined in
        the directions contained in this argument in order to converge. E.g.
        Grad_par requires refinement only in the y-direction.
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
    """


    error_list = []
    if full_test:
        print('testing',boutcore_operator)
    for n in ngrids:
        if full_test:
            print('n =',n)
        # set options
        boutcore.setOption('mxg', str(mxg), force=True)
        boutcore.setOption('myg', str(myg), force=True)
        # set up mesh input
        if 'x' in dimensions:
            boutcore.setOption('testmesh:nx', exprToStr(n+2*mxg), force=True)
        else:
            boutcore.setOption('testmesh:nx', exprToStr(default_n+2*mxg), force=True)
        if 'y' in dimensions:
            boutcore.setOption('testmesh:ny', exprToStr(n), force=True)
        else:
            boutcore.setOption('testmesh:ny', exprToStr(default_n), force=True)
        if 'z' in dimensions:
            boutcore.setOption('testmesh:nz', exprToStr(n), force=True)
        else:
            boutcore.setOption('testmesh:nz', exprToStr(default_n), force=True)
        boutcore.setOption('testmesh:dx', exprToStr(metric.psiwidth*metric.scalex/n), force=True)
        boutcore.setOption('testmesh:dy', exprToStr(2.*pi*metric.scaley/n), force=True)
        boutcore.setOption('testmesh:dz', exprToStr(2.*pi/n), force=True)
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

        # create new Mesh object
        mesh = boutcore.Mesh(section='testmesh')

        # calculate result of differential operator using BOUT++ implementation
        # make this_testfunc periodic in radial-poloidal-toroidal coordinates
        # (no branch cuts)
        this_testfunc = testfunc.subs(metric.z, metric.z + metric.zShift) # ensure we don't change global testfunc

        bout_input = boutcore.create3D(exprToStr(this_testfunc), mesh)

        mesh.communicate(bout_input)
        if method is None:
            bout_result = boutcore_operator(bout_input)
        else:
            bout_result = boutcore_operator(bout_input, method=method)

        # Calculate result of differential operator symbolically, then convert
        # to boutcore.Field3D/Field2D.
        # Only want analytic derivatives to include shift for y-derivatives,
        # not x-derivatives. Note this is taken care of by mms_alternate
        analytic_func = symbolic_operator(testfunc) # take derivatives

        # shift analytic_func to pass to BOUT++, since create3D will shift from
        # field-aligned to orthogonal
        analytic_func = analytic_func.subs(metric.z, metric.z + metric.zShift)
        analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh)

        # calculate max error
        error = bout_result - analytic_result # as Field3D/Field2D
        error = error.get()[mxg:-mxg, myg:-myg] # numpy array, without guard cells
        error_list.append(numpy.max(numpy.abs(error))) # max error

    logerrors = numpy.log(error_list[-2]/error_list[-1])
    logspacing = numpy.log(ngrids[-1]/ngrids[-2])
    convergence = logerrors/logspacing

    if order-.1 < convergence < order+.2:
        return ['pass']
    else:
        print(error_list)
        if plot_error:
            try:
                from matplotlib import pyplot
                pyplot.loglog(1./ngrids, error_list)
                pyplot.show()
                from boututils.showdata import showdata
                showdata([error, bout_result.get()[mxg:-mxg, myg:-myg], analytic_result.get()[mxg:-mxg, myg:-myg]], titles=["error", "BOUT++ result", "analytic result"])
            except:
                pass
        return [str(boutcore_operator)+' is not working for '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.']

def test_operator2(dimensions, ngrids, testfunc1, testfunc2, boutcore_operator, symbolic_operator, order, method=None):
    """
    Test a set of parameters

    Parameters
    ----------

    dimensions : str
        String containing any of 'x', 'y' and 'z': the grid must be refined in
        the directions contained in this argument in order to converge. E.g.
        Grad_par requires refinement only in the y-direction.
    ngrids : numpy.array(int)
        Array of grid sizes to use. All directions being refined are given the
        same size.
    testfunc1, testfunc2 : sympy expression
        The inputs that will be given to the function being tested
    boutcore_operator : function
        function from boutcore to be tested
    boutcore_operator : function
        function from boutcore to be tested
    symbolic_operator : function
        function using sympy to do the symbolic equivalent of boutcore_operator
    order : int
        expected order of convergence of boutcore_operator
    method : str
        DIFF_METHOD to pass to boutcore_operator. If method=None then method
        argument will not be passed, so the default will be used.
    """


    error_list = []
    if full_test:
        print('testing', boutcore_operator)
    for n in ngrids:
        if full_test:
            print('n =',n)
        # set options
        boutcore.setOption('mxg', str(mxg), force=True)
        boutcore.setOption('myg', str(myg), force=True)
        # set up mesh input
        if 'x' in dimensions:
            boutcore.setOption('testmesh:nx', exprToStr(n+2*mxg), force=True)
        else:
            boutcore.setOption('testmesh:nx', exprToStr(default_n+2*mxg), force=True)
        if 'y' in dimensions:
            boutcore.setOption('testmesh:ny', exprToStr(n), force=True)
        else:
            boutcore.setOption('testmesh:ny', exprToStr(default_n), force=True)
        if 'z' in dimensions:
            boutcore.setOption('testmesh:nz', exprToStr(n), force=True)
        else:
            boutcore.setOption('testmesh:nz', exprToStr(default_n), force=True)
        boutcore.setOption('testmesh:dx', exprToStr(metric.psiwidth*metric.scalex/n), force=True)
        boutcore.setOption('testmesh:dy', exprToStr(2.*pi*metric.scaley/n), force=True)
        boutcore.setOption('testmesh:dz', exprToStr(2.*pi/n), force=True)
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

        # create new Mesh object
        mesh = boutcore.Mesh(section='testmesh')

        # calculate result of differential operator using BOUT++ implementation
        # make this_testfunc1 periodic in radial-poloidal-toroidal coordinates
        # (no branch cuts)
        this_testfunc1 = testfunc1.subs(metric.z, metric.z + metric.zShift)

        bout_input1 = boutcore.create3D(exprToStr(this_testfunc1), mesh)

        # make this_testfunc2 periodic in radial-poloidal-toroidal coordinates
        # (no branch cuts)
        this_testfunc2 = testfunc2.subs(metric.z, metric.z + metric.zShift)

        bout_input2 = boutcore.create3D(exprToStr(this_testfunc2), mesh)

        mesh.communicate(bout_input1, bout_input2)
        if method is None:
            bout_result = boutcore_operator(bout_input1, bout_input2)
        else:
            bout_result = boutcore_operator(bout_input1, bout_input2, method=method)

        # Calculate result of differential operator symbolically, then convert
        # to boutcore.Field3D/Field2D.
        # Only want analytic derivatives to include shift for y-derivatives,
        # not x-derivatives. Note this is taken care of by mms_alternate
        analytic_func = symbolic_operator(testfunc1, testfunc2) # take derivatives

        # shift analytic_func to pass to BOUT++, since create3D will shift from
        # field-aligned to orthogonal
        analytic_func = analytic_func.subs(metric.z, metric.z + metric.zShift)
        analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh)

        # calculate max error
        error = bout_result - analytic_result # as Field3D
        error = error.get()[mxg:-mxg, myg:-myg] # numpy array, without guard cells
        error_list.append(numpy.max(numpy.abs(error))) # max error

    logerrors = numpy.log(error_list[-2]/error_list[-1])
    logspacing = numpy.log(ngrids[-1]/ngrids[-2])
    convergence = logerrors/logspacing

    if order-.1 < convergence < order+.2:
        return ['pass']
    else:
        if plot_error:
            try:
                from matplotlib import pyplot
                pyplot.loglog(1./ngrids, error_list)
                pyplot.show()
                from boututils.showdata import showdata
                showdata([error, bout_result.get()[mxg:-mxg, myg:-myg], analytic_result.get()[mxg:-mxg, myg:-myg]], titles=["error", "BOUT++ result", "analytic result"])
            except:
                pass
        return [str(boutcore_operator)+' is not working for '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.']

min_exponent = 6
max_exponent = 7
ngrids = numpy.logspace(min_exponent, max_exponent, num=max_exponent-min_exponent+1, base=2).astype('int')
default_n = 4
mxg = 2
myg = 2
testfunc = cos(2*pi*metric.x+metric.y+metric.z)
testfunc2 = sin(4*pi*metric.x+2*metric.y+2*metric.z)+cos(2*pi*metric.x-metric.z)
order = 2
plot_error = False
test_deriv_ops = full_test
tests_3d = full_test

boutcore.init('-q -q -q -q')

results = []

# single-argument operators
results += test_operator('yz', ngrids, testfunc, boutcore.Grad_par, Grad_par, order)
results += test_operator('yz', ngrids, testfunc, boutcore.Div_par, Div_par, order)
results += test_operator('yz', ngrids, testfunc, boutcore.Grad2_par2, Grad2_par2, order)
if tests_3d:
    results += test_operator('xyz', ngrids, testfunc, boutcore.Laplace, Laplace, order)
results += test_operator('yz', ngrids, testfunc, boutcore.Laplace_par, Laplace_par, order)
# note Laplace_perp uses Laplace, so needs y-dimension refinement to converge
if tests_3d:
    results += test_operator('xyz', ngrids, testfunc, boutcore.Laplace_perp, Laplace_perp, order)
# Delp2 uses the global mesh, which we can't reset, so can't test here
#results += test_operator('xz', ngrids, testfunc, boutcore.Delp2, Delp2, order)

# two-argument operators
results += test_operator2('yz', ngrids, testfunc, testfunc2, boutcore.Vpar_Grad_par, Vpar_Grad_par, 1, method='U1') # can only use U1 upwinding with ShiftedMetric (at least until 2 yup/ydown fields are implemented)
results += test_operator2('yz', ngrids, testfunc, testfunc2, boutcore.Div_par_K_Grad_par, Div_par_K_Grad_par, order)
results += test_operator2('yz', ngrids, testfunc, testfunc2, boutcore.Div_par_flux, lambda v,f: Div_par(v*f), 1)
results += test_operator2('yz', ngrids, testfunc, testfunc2, boutcore.Div_par_flux, lambda v,f: Div_par(v*f), order, method='C2')
# note bracket(Field2D, Field2D) is exactly zero, so doesn't make sense to MMS test
# Note BRACKET_STD version of bracket includes parallel derivatives, so needs
# y-dimension refinement to converge.
if tests_3d:
    results += test_operator2('xyz', ngrids, testfunc, testfunc2, boutcore.bracket, lambda a,b: b0xGrad_dot_Grad(a,b)/metric.B, order, method='BRACKET_STD')

if test_deriv_ops:
    # test derivative operators
    results += test_operator('yz', ngrids, testfunc, boutcore.DDY, DDY, order)
    results += test_operator('yz', ngrids, testfunc, boutcore.D2DY2, D2DY2, order)
    # D4DY4 requires 5-point stencil, so can't test with ShiftedMetric
    #results += test_operator('yz', ngrids, testfunc, boutcore.D4DY4, D4DY4, order)
    results += test_operator2('yz', ngrids, testfunc, testfunc2, boutcore.VDDY, lambda v,f: v*DDY(f), 2)
    results += test_operator2('yz', ngrids, testfunc, testfunc2, boutcore.FDDY, lambda v,f: DDY(v*f), 1)
    if tests_3d:
        results += test_operator('xyz', ngrids, testfunc, boutcore.D2DXDY, D2DXDY, order)
        results += test_operator('xyz', ngrids, testfunc, boutcore.D2DYDX, D2DYDX, order)
    results += test_operator('yz', ngrids, testfunc, boutcore.D2DYDZ, D2DYDZ, order)
    results += test_operator('yz', ngrids, testfunc, boutcore.D2DZDY, D2DZDY, order)

# check results of tests
fail = False
for result in results:
    if result is not 'pass':
        print(result)
        fail = True
if fail:
    exit(1)
else:
    print('pass')
    exit(0)
