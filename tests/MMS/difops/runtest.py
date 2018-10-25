#!/usr/bin/env python3

# MMS test for differential operators (that use the metric)

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
parser.add_argument('--operator', default=None)
args = parser.parse_args()
full_test = not args.short

min_exponent = 6
max_exponent = 7
ngrids = numpy.logspace(min_exponent, max_exponent, num=max_exponent-min_exponent+1, base=2).astype(int)
ngrids_lowres = numpy.logspace(6, 7, num=2, base=2).astype(int)
default_n = 1
mxg = 2
myg = 2
testfunc = cos(2*pi*metric.x+metric.y+metric.z)
testfunc2 = sin(4*pi*metric.x+2*metric.y+2*metric.z)+cos(2*pi*metric.x-metric.z)
accuracy_order = 2
plot_error = False
test_deriv_ops = full_test
tests_3d = full_test
test_throw = full_test

# geometry for simple circular tokamak
# Domain in SOL to test y-boundary conditions
tokamak = SimpleTokamak(psiN0=1.1)

# rescale x and y coordinates so dx and dy are not constants
tokamak.set_scalex(1 + .1*sin(2*pi*metric.x+metric.y))
tokamak.set_scaley(1 + .1*sin(2*pi*metric.x-metric.y))
# re-calculate metric terms
tokamak.metric()

# single-argument operators
all_types = ('2D', '3D') # eventually, should be able to use this when Field2D operators support staggering properly
type_3d = ('3D',)
type_2d = ('2D',)

# dict of tuples contaning: (<directions to stagger>, <directions to refine grid>, <array of grid sizes>, <input function>,
#                            <boutcore operator>, <mms_alternate symbolic operator>, <order of (highest) derivative>,
#                            <expected order of accuracy>, <tuple of field types to test>, <method [optional])
single_argument_operators = {
    'Grad_par': ('y', 'y', ngrids, testfunc, boutcore.Grad_par, Grad_par, 1, accuracy_order, type_3d), # staggering in y-direction allowed
    'Grad_par_2d': ('', 'y', ngrids, testfunc, boutcore.Grad_par, Grad_par, 1, accuracy_order, type_2d), # no staggering allowed
    'Div_par': ('y', 'y', ngrids, testfunc, boutcore.Div_par, Div_par, 1, accuracy_order, type_3d), # staggering in y-direction allowed
    'Div_par_2d': ('', 'y', ngrids, testfunc, boutcore.Div_par, Div_par, 1, accuracy_order, type_2d), # no staggering allowed
    'Grad2_par2': ('y', 'y', ngrids, testfunc, boutcore.Grad2_par2, Grad2_par2, 2, accuracy_order, type_3d), # staggering in y-direction allowed
    'Grad2_par2_2d': ('', 'y', ngrids, testfunc, boutcore.Grad2_par2, Grad2_par2, 2, accuracy_order, type_2d), # no staggering allowed
    'Laplace': ('', 'xyz', ngrids_lowres, testfunc, boutcore.Laplace, Laplace, 2, accuracy_order, all_types), # no staggering allowed
    'Laplace_par': ('y', 'y', ngrids, testfunc, boutcore.Laplace_par, Laplace_par, 2, accuracy_order, type_3d), # staggering in y-direction allowed
    'Laplace_par_2d': ('', 'y', ngrids, testfunc, boutcore.Laplace_par, Laplace_par, 2, accuracy_order, type_2d), # no staggering allowed
    # note Laplace_perp uses Laplace, so needs y-dimension refinement to converge
    'Laplace_perp': ('', 'xyz', ngrids_lowres, testfunc, boutcore.Laplace_perp, Laplace_perp, 2, accuracy_order, all_types), # no staggering allowed
    # Delp2 uses the global mesh, which we can't reset, so can't test here
    #'Delp2': ('x', 'xz', ngrids, testfunc, boutcore.Delp2, Delp2, accuracy_order),
}

# two-argument operators
all_types2 = [('2D', '2D'), ('3D', '3D')] # expand this to include mixed 2D/3D types at some point
types_2d = [('2D', '2D')]
types_3d = [('3D', '3D')]

# dict of tuples contaning: (<directions to stagger>, <directions to refine grid>, <array of grid sizes>, <first input function>,
#                            <second input function>, <boutcore operator>, <mms_alternate symbolic operator>, <order of (highest) derivative>,
#                            <expected order of accuracy>, <tuple of field types to test>, <method [optional]>)
two_argument_operators = {
    'Vpar_Grad_par': ('y', 'y', ngrids, testfunc, testfunc2, boutcore.Vpar_Grad_par, Vpar_Grad_par, 1, accuracy_order, types_3d), # some staggering in y-direction allowed
    'Vpar_Grad_par_2d': ('y', 'y', ngrids, testfunc, testfunc2, boutcore.Vpar_Grad_par, Vpar_Grad_par, 1, accuracy_order, types_2d), # no staggering allowed
    'Div_par_K_Grad_par': ('yy', 'y', ngrids, testfunc, testfunc2, boutcore.Div_par_K_Grad_par, Div_par_K_Grad_par, 1, accuracy_order, types_3d), # any staggering in y-direction allowed
    'Div_par_K_Grad_par_2d': ('', 'y', ngrids, testfunc, testfunc2, boutcore.Div_par_K_Grad_par, Div_par_K_Grad_par, 1, accuracy_order, types_2d), # no staggering allowed
    'Div_par_flux': ('y', 'y', ngrids, testfunc, testfunc2, boutcore.Div_par_flux, lambda v,f: Div_par(v*f), 1, 1, types_3d), # some staggering in y-direction allowed
    'Div_par_flux_2d': ('', 'y', ngrids, testfunc, testfunc2, boutcore.Div_par_flux, lambda v,f: Div_par(v*f), 1, 1, types_2d), # no staggering allowed
    'Div_par_flux_C2': ('y', 'y', ngrids, testfunc, testfunc2, boutcore.Div_par_flux, lambda v,f: Div_par(v*f), 1, accuracy_order, types_3d, 'C2'), # some staggering in y-direction allowed
    'Div_par_flux_C2_2d': ('', 'y', ngrids, testfunc, testfunc2, boutcore.Div_par_flux, lambda v,f: Div_par(v*f), 1, accuracy_order, types_2d, 'C2'), # no staggering allowed
    # note bracket(Field2D, Field2D) is exactly zero, so doesn't make sense to MMS test
    'bracket_ARAKAWA': ('', 'xz', ngrids, testfunc, testfunc2, boutcore.bracket, bracket, 1, accuracy_order, types_3d, 'BRACKET_ARAKAWA'), # no staggering allowed
    # Note BRACKET_STD version of bracket includes parallel derivatives, so needs
    # y-dimension refinement to converge.
    # Also it converges faster than 2nd order at 64->128, but approaches closer
    # when resolution is increased. Allow test to pass anyway by increasing
    # expected order
    'bracket_STD': ('', 'xyz', ngrids_lowres, testfunc, testfunc2, boutcore.bracket, lambda a,b: b0xGrad_dot_Grad(a,b)/metric.B, 1, 2.6, types_3d, 'BRACKET_STD'), # no staggering allowed
}

# test derivative operators
# dict of tuples contaning: (<directions to stagger>, <directions to refine grid>, <array of grid sizes>, <input function>,
#                            <boutcore operator>, <mms_alternate symbolic operator>, <order of (highest) derivative>,
#                            <expected order of accuracy>, <tuple of field types to test>, <method [optional])
single_argument_derivative_operators = {
    'DDX': ('x', 'x', ngrids, testfunc, boutcore.DDX, DDX, 1, accuracy_order, type_3d),
    'DDX_2d': ('', 'x', ngrids, testfunc, boutcore.DDX, DDX, 1, accuracy_order, type_2d),
    'DDY': ('y', 'y', ngrids, testfunc, boutcore.DDY, DDY, 1, accuracy_order, type_3d),
    'DDY_2d': ('', 'y', ngrids, testfunc, boutcore.DDY, DDY, 1, accuracy_order, type_2d),
    'DDZ': ('', 'z', ngrids, testfunc, boutcore.DDZ, DDZ, 1, accuracy_order, type_3d, 'C2'),
    'D2DX2': ('x', 'x', ngrids, testfunc, boutcore.D2DX2, D2DX2, 2, accuracy_order, type_3d),
    'D2DX2_2d': ('', 'x', ngrids, testfunc, boutcore.D2DX2, D2DX2, 2, accuracy_order, type_2d),
    'D2DY2': ('y', 'y', ngrids, testfunc, boutcore.D2DY2, D2DY2, 2, accuracy_order, type_3d),
    'D2DY2_2d': ('', 'y', ngrids, testfunc, boutcore.D2DY2, D2DY2, 2, accuracy_order, type_2d),
    'D2DZ2': ('', 'z', ngrids, testfunc, boutcore.D2DZ2, D2DZ2, 2, accuracy_order, type_3d, 'C2'),
    'D4DX4': ('', 'x', ngrids, testfunc, boutcore.D4DX4, D4DX4, 4, accuracy_order, type_3d),
    'D4DX4_2d': ('', 'x', ngrids, testfunc, boutcore.D4DX4, D4DX4, 4, accuracy_order, type_2d),
    'D4DY4': ('', 'y', ngrids, testfunc, boutcore.D4DY4, D4DY4, 4, accuracy_order, type_3d),
    'D4DY4_2d': ('', 'y', ngrids, testfunc, boutcore.D4DY4, D4DY4, 4, accuracy_order, type_2d),
    'D4DZ4': ('', 'z', ngrids, testfunc, boutcore.D4DZ4, D4DZ4, 4, accuracy_order, type_3d), # D4DZ4 is hard coded to use DIFF_C2
    'D2DXDY': ('y', 'xy', ngrids, testfunc, boutcore.D2DXDY, D2DXDY, 2, accuracy_order, type_3d),
    'D2DXDY_2d': ('', 'xy', ngrids, testfunc, boutcore.D2DXDY, D2DXDY, 2, accuracy_order, type_2d),
    'D2DYDX': ('x', 'xy', ngrids, testfunc, boutcore.D2DYDX, D2DYDX, 2, accuracy_order, type_3d),
    'D2DYDX_2d': ('', 'xy', ngrids, testfunc, boutcore.D2DYDX, D2DYDX, 2, accuracy_order, type_2d),
    'D2DXDZ': ('', 'xz', ngrids, testfunc, boutcore.D2DXDZ, D2DXDZ, 2, accuracy_order, type_3d),
    'D2DZDX': ('x', 'xz', ngrids, testfunc, boutcore.D2DZDX, D2DZDX, 2, accuracy_order, type_3d),
    'D2DYDZ': ('', 'yz', ngrids, testfunc, boutcore.D2DYDZ, D2DYDZ, 2, accuracy_order, type_3d),
    'D2DZDY': ('y', 'yz', ngrids, testfunc, boutcore.D2DZDY, D2DZDY, 2, accuracy_order, type_3d),
}

# dict of tuples contaning: (<directions to stagger>, <directions to refine grid>, <array of grid sizes>, <first input function>,
#                            <second input function>, <boutcore operator>, <mms_alternate symbolic operator>, <order of (highest) derivative>,
#                            <expected order of accuracy>, <tuple of field types to test>, <method [optional]>)
two_argument_derivative_operators = {
    'VDDX': ('x', 'x', ngrids, testfunc, testfunc2, boutcore.VDDX, lambda v,f: v*DDX(f), 1, accuracy_order, types_3d),
    'VDDX_2d': ('', 'x', ngrids, testfunc, testfunc2, boutcore.VDDX, lambda v,f: v*DDX(f), 1, accuracy_order, types_2d),
    'VDDY': ('y', 'y', ngrids, testfunc, testfunc2, boutcore.VDDY, lambda v,f: v*DDY(f), 1, accuracy_order, types_3d),
    'VDDY_2d': ('y', 'y', ngrids, testfunc, testfunc2, boutcore.VDDY, lambda v,f: v*DDY(f), 1, accuracy_order, types_2d),
    'VDDZ': ('z', 'z', ngrids, testfunc, testfunc2, boutcore.VDDZ, lambda v,f: v*DDZ(f), 1, accuracy_order, types_3d, 'C2'),
    'FDDX': ('x', 'x', ngrids, testfunc, testfunc2, boutcore.FDDX, lambda v,f: DDX(v*f), 1, 1, types_3d),
    'FDDX_2d': ('', 'x', ngrids, testfunc, testfunc2, boutcore.FDDX, lambda v,f: DDX(v*f), 1, 1, types_2d),
    'FDDY': ('y', 'y', ngrids, testfunc, testfunc2, boutcore.FDDY, lambda v,f: DDY(v*f), 1, 1, types_3d),
    'FDDY_2d': ('', 'y', ngrids, testfunc, testfunc2, boutcore.FDDY, lambda v,f: DDY(v*f), 1, 1, types_2d),
    'FDDZ': ('z', 'z', ngrids, testfunc, testfunc2, boutcore.FDDZ, lambda v,f: DDZ(v*f), 1, 1, types_3d),
}

if full_test:
    boundary_conditions = ["", "dirichlet_o", "neumann_o", "free_o"]
else:
    boundary_conditions = [""]

def test_operator(ngrids, testfunc, dimensions, boutcore_operator, symbolic_operator, accuracy_order, ftype, method, stagger, boundary_condition):

    testfunc = copy(testfunc) # ensure we don't change global testfunc
    error_list = []
    if full_test:
        print('testing',boutcore_operator, ftype, stagger, boundary_condition)

    for n in ngrids:
        if full_test:
            print('n =',n)
        # set options
        # set up mesh input
        if 'x' in dimensions:
            this_mxg = mxg
            nx = n+2*this_mxg
        else:
            this_mxg = 1
            nx = default_n+2
        boutcore.setOption('mxg', str(this_mxg), force=True)
        boutcore.setOption('testmesh:nx', exprToStr(nx), force=True)
        if 'y' in dimensions:
            this_myg = myg
            ny = n
        else:
            this_myg = 1
            ny = default_n
        boutcore.setOption('myg', str(this_myg), force=True)
        boutcore.setOption('testmesh:ny', exprToStr(ny), force=True)
        if 'z' in dimensions:
            nz = n
        else:
            nz = default_n
        boutcore.setOption('testmesh:nz', exprToStr(nz), force=True)

        boutcore.setOption('testmesh:dx', exprToStr(metric.psiwidth*metric.scalex/nx), force=True)
        boutcore.setOption('testmesh:dy', exprToStr(2.*pi*metric.scaley/ny), force=True)
        boutcore.setOption('testmesh:dz', exprToStr(2.*pi/nz), force=True)
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
        tokamak.setNs(nx=nx, ny=ny, nz=nz)
        boutcore.setOption('testmesh:ixseps1', exprToStr(tokamak.ixseps1), force=True)
        boutcore.setOption('testmesh:ixseps2', exprToStr(tokamak.ixseps2), force=True)
        if ny>1:
            boutcore.setOption('testmesh:jyseps1_1', exprToStr(tokamak.jyseps1_1), force=True)
            boutcore.setOption('testmesh:jyseps1_2', exprToStr(tokamak.jyseps1_2), force=True)
            boutcore.setOption('testmesh:jyseps2_1', exprToStr(tokamak.jyseps2_1), force=True)
            boutcore.setOption('testmesh:jyseps2_2', exprToStr(tokamak.jyseps2_2), force=True)
        else:
            boutcore.setOption('testmesh:jyseps1_1', exprToStr(-1), force=True)
            boutcore.setOption('testmesh:jyseps1_2', exprToStr(0), force=True)
            boutcore.setOption('testmesh:jyseps2_1', exprToStr(0), force=True)
            boutcore.setOption('testmesh:jyseps2_2', exprToStr(2), force=True)
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
        if ftype == '2D':
            # cannot have z-dependence
            testfunc = testfunc.subs(metric.z, 0)
            bout_input = boutcore.create2D(exprToStr(testfunc), mesh, outloc=inloc)
        elif ftype == '3D':
            bout_input = boutcore.create3D(exprToStr(testfunc), mesh, outloc=inloc)
        else:
            raise ValueError('Unexpected ftype argument '+str(ftype))
        if boundary_condition is "":
            pass
        elif "dirichlet" in boundary_condition:
            if nx>3:
                bout_input.applyBoundary("%s(%s)"%(boundary_condition, exprToStr(testfunc)), region="core")
                bout_input.applyBoundary("%s(%s)"%(boundary_condition, exprToStr(testfunc)), region="sol")
            if ny>1:
                bout_input.applyBoundary("%s(%s)"%(boundary_condition, exprToStr(testfunc)), region="lower_target")
                bout_input.applyBoundary("%s(%s)"%(boundary_condition, exprToStr(testfunc)), region="upper_target")
        elif "neumann" in boundary_condition:
            if nx>3:
                bout_input.applyBoundary(boundary="%s(%s)"%(boundary_condition, exprToStr(DDX(testfunc))), region="core")
                bout_input.applyBoundary(boundary="%s(%s)"%(boundary_condition, exprToStr(DDX(testfunc))), region="sol")
            if ny>1:
                bout_input.applyBoundary(boundary="%s(%s)"%(boundary_condition, exprToStr(DDY(testfunc))), region="lower_target")
                bout_input.applyBoundary(boundary="%s(%s)"%(boundary_condition, exprToStr(DDY(testfunc))), region="upper_target")
        else:
            # free boundary, does not need argument
            if nx>3:
                bout_input.applyBoundary(boundary_condition, region="core")
                bout_input.applyBoundary(boundary_condition, region="sol")
            if ny>1:
                bout_input.applyBoundary(boundary_condition, region="lower_target")
                bout_input.applyBoundary(boundary_condition, region="upper_target")
        if method is None:
            bout_result = boutcore_operator(bout_input, outloc=outloc)
        else:
            bout_result = boutcore_operator(bout_input, outloc=outloc, method=method)

        # calculate result of differential operator symbolically, then convert to boutcore.Field3D/Field2D
        analytic_func = symbolic_operator(testfunc)
        if ftype == '2D':
            analytic_result = boutcore.create2D(exprToStr(analytic_func), mesh, outloc=outloc)
        elif ftype == '3D':
            analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh, outloc=outloc)
        else:
            raise ValueError('Unexpected ftype argument '+str(ftype))

        # calculate max error
        error = bout_result - analytic_result # as Field3D/Field2D
        error = error.get()[this_mxg:-this_mxg, this_myg:-this_myg] # numpy array, without guard cells
        error_list.append(numpy.max(numpy.abs(error))) # max error

    logerrors = numpy.log(error_list[-2]/error_list[-1])
    logspacing = numpy.log(ngrids[-1]/ngrids[-2])
    convergence = logerrors/logspacing

    if accuracy_order-.1 < convergence < accuracy_order+.2:
        return ['pass']
    else:
        if plot_error:
            from matplotlib import pyplot
            pyplot.loglog(1./ngrids, error_list)
            pyplot.show()
            from boututils.showdata import showdata
            showdata([numpy.atleast_2d(numpy.squeeze(error)),
                numpy.atleast_2d(numpy.squeeze(bout_result.get()[this_mxg:-this_mxg, this_myg:-this_myg])),
                numpy.atleast_2d(numpy.squeeze(analytic_result.get()[this_mxg:-this_mxg, this_myg:-this_myg]))])
        print(str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' '+str(ftype)+' '+str(method)+' '+str(boundary_condition)+'. Expected '+str(accuracy_order)+', got '+str(convergence)+'.')
        return [str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' '+str(ftype)+' '+str(method)+' '+str(boundary_condition)+'. Expected '+str(accuracy_order)+', got '+str(convergence)+'.']

def cycle_staggering(stagger_directions, base_dimensions, ngrids, testfunc, boutcore_operator, symbolic_operator, deriv_order, accuracy_order, types, method=None):
    """
    Loop over different parameters, calling test_operator for each

    Parameters
    ----------

    stagger_directions : str
        Directions in which this operator can be staggered. String containing
        'x', 'y' or 'z' will test both permutations of staggering between
        CENTRE and the corresponding direction.
    base_dimensions : str
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
    types : list of str
        list of strings giving the type of the argument to the operators
        respectively. '2D' for Field2D and '3D' for Field3D.
    """

    # all derivatives at same inloc/outloc should work
    dimensions_staggers = [(base_dimensions, None), # no staggering
                           (base_dimensions+'', ('CENTRE', 'CENTRE')), # include all-centred, but with StaggerGrids=true
                           (base_dimensions+'x', ('XLOW', 'XLOW')),
                           (base_dimensions+'y', ('YLOW', 'YLOW')),
                           (base_dimensions+'z', ('ZLOW', 'ZLOW'))]
    if 'x' in stagger_directions:
        dimensions_staggers += [(base_dimensions+'x', ('CENTRE', 'XLOW')),
                                (base_dimensions+'x', ('XLOW', 'CENTRE'))]
    if 'y' in stagger_directions:
        dimensions_staggers += [(base_dimensions+'y', ('CENTRE', 'YLOW')),
                                (base_dimensions+'y', ('YLOW', 'CENTRE'))]
    if 'z' in stagger_directions:
        dimensions_staggers += [(base_dimensions+'z', ('CENTRE', 'ZLOW')),
                                (base_dimensions+'z', ('ZLOW', 'CENTRE'))]
    result = []
    for boundary_condition in boundary_conditions:
        if 'dirichlet' in boundary_condition or 'free' in boundary_condition:
            if deriv_order == 1:
                boundary_condition += '3'
            elif deriv_order == 2:
                boundary_condition += '4'
            else:
                boundary_condition += '5' # dirichlet_o5/free_o5 is highest order available
        elif 'neumann' in boundary_condition:
            if deriv_order == 1:
                boundary_condition += '2'
            elif deriv_order == 2:
                boundary_condition += '4' # neumann_o3 is not implemented
            else:
                boundary_condition += '4' # neumann_o4 is highest order available
        for ftype in types:
            for dimensions, stagger in dimensions_staggers:
                if tests_3d or not ('x' in dimensions and 'y' in dimensions and 'z' in dimensions):
                    result += test_operator(ngrids, testfunc, dimensions, boutcore_operator, symbolic_operator, accuracy_order, ftype, method, stagger, boundary_condition)

            if test_throw:
                # check that unsupported combinations of locations throw an exception

                locations = ['CENTRE', 'XLOW', 'YLOW', 'ZLOW']
                # first make a list of all permutations
                fail_staggers = [(x,y) for x in locations for y in locations]
                # now remove the ones have already tested
                for dimensions, stagger in dimensions_staggers:
                    if stagger is not None:
                        index = fail_staggers.index(stagger)
                        del fail_staggers[index]
                if boundary_condition is "":
                    # Doesn't add anything to check for all boundary conditions
                    for stagger in fail_staggers:
                        # check that an exception is throw for combinations of directions that we expect to fail
                        try:
                            test = test_operator(numpy.array([16, 16]), testfunc, 'xyz', boutcore_operator, symbolic_operator, accuracy_order, ftype, method, stagger, boundary_condition)
                        except RuntimeError:
                            result += ['pass']
                        else:
                            result += ['Expected '+str(boutcore_operator)+' to throw for '+stagger[0]+'->'+stagger[1]+' '+' '+str(ftype)+' '+str(method)+' but it did not.']

    return result

def test_operator2(ngrids, testfunc1, testfunc2, dimensions, boutcore_operator, symbolic_operator, accuracy_order, ftypes, method, stagger):

    testfunc1 = copy(testfunc1) # ensure we don't change global testfunc1
    testfunc2 = copy(testfunc2) # ensure we don't change global testfunc2
    error_list = []
    if full_test:
        print('testing', boutcore_operator, ftypes, stagger)
    for n in ngrids:
        if full_test:
            print('n =',n)
        # set options
        # set up mesh input
        if 'x' in dimensions:
            this_mxg = mxg
            nx = n+2*this_mxg
        else:
            this_mxg = 1
            nx = default_n+2*this_mxg
        boutcore.setOption('mxg', str(this_mxg), force=True)
        boutcore.setOption('testmesh:nx', exprToStr(nx), force=True)
        if 'y' in dimensions:
            this_myg = myg
            ny = n
        else:
            this_myg = 1
            ny = default_n
        boutcore.setOption('myg', str(this_myg), force=True)
        boutcore.setOption('testmesh:ny', exprToStr(ny), force=True)
        if 'z' in dimensions:
            nz = n
        else:
            nz = default_n
        boutcore.setOption('testmesh:nz', exprToStr(nz), force=True)

        boutcore.setOption('testmesh:dx', exprToStr(metric.psiwidth*metric.scalex/nx), force=True)
        boutcore.setOption('testmesh:dy', exprToStr(2.*pi*metric.scaley/ny), force=True)
        boutcore.setOption('testmesh:dz', exprToStr(2.*pi/nz), force=True)
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
        tokamak.setNs(nx=nx, ny=ny, nz=nz)
        boutcore.setOption('testmesh:ixseps1', exprToStr(tokamak.ixseps1), force=True)
        boutcore.setOption('testmesh:ixseps2', exprToStr(tokamak.ixseps2), force=True)
        if ny>1:
            boutcore.setOption('testmesh:jyseps1_1', exprToStr(tokamak.jyseps1_1), force=True)
            boutcore.setOption('testmesh:jyseps1_2', exprToStr(tokamak.jyseps1_2), force=True)
            boutcore.setOption('testmesh:jyseps2_1', exprToStr(tokamak.jyseps2_1), force=True)
            boutcore.setOption('testmesh:jyseps2_2', exprToStr(tokamak.jyseps2_2), force=True)
        else:
            boutcore.setOption('testmesh:jyseps1_1', exprToStr(-1), force=True)
            boutcore.setOption('testmesh:jyseps1_2', exprToStr(0), force=True)
            boutcore.setOption('testmesh:jyseps2_1', exprToStr(0), force=True)
            boutcore.setOption('testmesh:jyseps2_2', exprToStr(2), force=True)
        if stagger is None:
            boutcore.setOption('testmesh:staggergrids', str('false'), force=True)
        else:
            boutcore.setOption('testmesh:staggergrids', str('true'), force=True)

        # create new Mesh object
        mesh = boutcore.Mesh(section='testmesh')

        if stagger is None:
            vloc = 'CENTRE'
            inloc = 'CENTRE'
            outloc = 'CENTRE'
        else:
            vloc = stagger[0]
            inloc = stagger[1]
            outloc = stagger[2]

        # calculate result of differential operator using BOUT++ implementation
        if ftypes[0] == '2D':
            # cannot have z-dependence
            testfunc1 = testfunc1.subs(metric.z, 0)
            bout_input1 = boutcore.create2D(exprToStr(testfunc1), mesh, outloc=vloc)
        elif ftypes[0] == '3D':
            bout_input1 = boutcore.create3D(exprToStr(testfunc1), mesh, outloc=vloc)
        else:
            raise ValueError('Unexpected ftype argument '+str(ftypes[0]))
        if ftypes[1] == '2D':
            # cannot have z-dependence
            testfunc2 = testfunc2.subs(metric.z, 0)
            bout_input2 = boutcore.create2D(exprToStr(testfunc2), mesh, outloc=inloc)
        elif ftypes[1] == '3D':
            bout_input2 = boutcore.create3D(exprToStr(testfunc2), mesh, outloc=inloc)
        else:
            raise ValueError('Unexpected ftype argument '+str(ftypes[1]))
        if method is None:
            bout_result = boutcore_operator(bout_input1, bout_input2, outloc=outloc)
        else:
            bout_result = boutcore_operator(bout_input1, bout_input2, outloc=outloc, method=method)

        # calculate result of differential operator symbolically, then convert to boutcore.Field3D/Field2D
        analytic_func = symbolic_operator(testfunc1, testfunc2)
        if ftypes[0] == '2D' and ftypes[1] == '2D':
            analytic_result = boutcore.create2D(exprToStr(analytic_func), mesh, outloc=outloc)
        else:
            analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh, outloc=outloc)

        # calculate max error
        error = bout_result - analytic_result # as Field3D
        error = error.get()[this_mxg:-this_mxg, this_myg:-this_myg] # numpy array, without guard cells
        error_list.append(numpy.max(numpy.abs(error))) # max error

    logerrors = numpy.log(error_list[-2]/error_list[-1])
    logspacing = numpy.log(ngrids[-1]/ngrids[-2])
    convergence = logerrors/logspacing

    if accuracy_order-.1 < convergence < accuracy_order+.2:
        return ['pass']
    else:
        if plot_error:
            from matplotlib import pyplot
            pyplot.loglog(1./ngrids, error_list)
            pyplot.show()
            from boututils.showdata import showdata
            showdata([numpy.atleast_2d(numpy.squeeze(error)),
                numpy.atleast_2d(numpy.squeeze(bout_result.get()[this_mxg:-this_mxg, this_myg:-this_myg])),
                numpy.atleast_2d(numpy.squeeze(analytic_result.get()[this_mxg:-this_mxg, this_myg:-this_myg]))])
        print(str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' '+str(ftypes)+' '+str(method)+'. Expected '+str(accuracy_order)+', got '+str(convergence)+'.')
        return [str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' '+str(ftypes)+' '+str(method)+'. Expected '+str(accuracy_order)+', got '+str(convergence)+'.']

def cycle_staggering2(stagger_directions, base_dimensions, ngrids, testfunc1, testfunc2, boutcore_operator, symbolic_operator, deriv_order, accuracy_order, types, method=None):
    """
    Loop over different parameters, calling test_operator2 for each

    Parameters
    ----------

    stagger_directions : str
        Directions in which this operator can be staggered. String containing
        'xx', 'yy' or 'zz' will test all permutations of staggering between
        CENTRE and the corresponding direction. String containing 'x', 'y' or
        'z' will keep second argument and outloc at the same location, and
        stagger the first argument in the corresponding direction.
    base_dimensions : str
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
    symbolic_operator : function
        function using sympy to do the symbolic equivalent of boutcore_operator
    order : int
        expected order of convergence of boutcore_operator
    types : list of tuples of str
        list of pairs (type1, type2) giving the type of first and second
        arguments to the operators respectively. '2D' for Field2D and '3D' for
        Field3D.
    method : str
        DIFF_METHOD to pass to boutcore_operator. If method=None then method
        argument will not be passed, so the default will be used.
    """

    # all derivatives at same inloc/outloc should work
    dimensions_staggers = [(base_dimensions, None), # no staggering
                           (base_dimensions+'', ('CENTRE', 'CENTRE', 'CENTRE')), # include all-centred, but with StaggerGrids=true
                           (base_dimensions+'x', ('XLOW', 'XLOW', 'XLOW')),
                           (base_dimensions+'y', ('YLOW', 'YLOW', 'YLOW')),
                           (base_dimensions+'z', ('ZLOW', 'ZLOW', 'ZLOW'))]
    if 'xx' in stagger_directions:
        # for xx include all permutations of XLOW and CENTRE
        dimensions_staggers += [(base_dimensions+'x', ('CENTRE', 'CENTRE', 'XLOW')),
                                (base_dimensions+'x', ('CENTRE', 'XLOW', 'CENTRE')),
                                (base_dimensions+'x', ('XLOW', 'CENTRE', 'CENTRE')),
                                (base_dimensions+'x', ('XLOW', 'XLOW', 'CENTRE')),
                                (base_dimensions+'x', ('XLOW', 'CENTRE', 'XLOW')),
                                (base_dimensions+'x', ('CENTRE', 'XLOW', 'XLOW'))]
    elif 'x' in stagger_directions:
        dimensions_staggers += [(base_dimensions+'x', ('CENTRE', 'XLOW', 'XLOW')),
                                (base_dimensions+'x', ('XLOW', 'CENTRE', 'CENTRE'))]
    if 'yy' in stagger_directions:
        # for yy include all permutations of YLOW and CENTRE
        dimensions_staggers += [(base_dimensions+'y', ('CENTRE', 'CENTRE', 'YLOW')),
                                (base_dimensions+'y', ('CENTRE', 'YLOW', 'CENTRE')),
                                (base_dimensions+'y', ('YLOW', 'CENTRE', 'CENTRE')),
                                (base_dimensions+'y', ('YLOW', 'YLOW', 'CENTRE')),
                                (base_dimensions+'y', ('YLOW', 'CENTRE', 'YLOW')),
                                (base_dimensions+'y', ('CENTRE', 'YLOW', 'YLOW'))]
    elif 'y' in stagger_directions:
        dimensions_staggers += [(base_dimensions+'y', ('CENTRE', 'YLOW', 'YLOW')),
                                (base_dimensions+'y', ('YLOW', 'CENTRE', 'CENTRE'))]
    if 'zz' in stagger_directions:
        # for zz include all permutations of ZLOW and CENTRE
        dimensions_staggers += [(base_dimensions+'z', ('CENTRE', 'CENTRE', 'ZLOW')),
                                (base_dimensions+'z', ('CENTRE', 'ZLOW', 'CENTRE')),
                                (base_dimensions+'z', ('ZLOW', 'CENTRE', 'CENTRE')),
                                (base_dimensions+'z', ('ZLOW', 'ZLOW', 'CENTRE')),
                                (base_dimensions+'z', ('ZLOW', 'CENTRE', 'ZLOW')),
                                (base_dimensions+'z', ('CENTRE', 'ZLOW', 'ZLOW'))]
    elif 'z' in stagger_directions:
        dimensions_staggers += [(base_dimensions+'z', ('CENTRE', 'ZLOW', 'ZLOW')),
                                (base_dimensions+'z', ('ZLOW', 'CENTRE', 'CENTRE'))]
    result = []
    for ftypes in types:
        for dimensions, stagger in dimensions_staggers:
            if tests_3d or not ('x' in dimensions and 'y' in dimensions and 'z' in dimensions):
                result += test_operator2(ngrids, testfunc1, testfunc2, dimensions, boutcore_operator, symbolic_operator, accuracy_order, ftypes, method, stagger)

        if test_throw:
            # check that unsupported combinations of locations throw an exception

            locations = ['CENTRE', 'XLOW', 'YLOW', 'ZLOW']
            # first make a list of all permutations
            fail_staggers = [(x,y,z) for x in locations for y in locations for z in locations]
            # now remove the ones have already tested
            for dimensions, stagger in dimensions_staggers:
                if stagger is not None:
                    index = fail_staggers.index(stagger)
                    del fail_staggers[index]
            for stagger in fail_staggers:
                # check that an exception is throw for combinations of directions that we expect to fail
                try:
                    test = test_operator2(numpy.array([16, 16]), testfunc1, testfunc2, 'xyz', boutcore_operator, symbolic_operator, accuracy_order, ftypes, method, stagger)
                except RuntimeError:
                    result += ['pass']
                else:
                    result += ['Expected '+str(boutcore_operator)+' to throw for '+stagger[0]+','+stagger[1]+'->'+stagger[2]+' '+str(ftypes)+' '+str(method)+' but it did not.']

    return result

if test_throw:
    if boutcore.bout_CHECK < 1:
        print('Warning: CHECK='+str(boutcore.bout_CHECK)+' so exceptions will '
            'not be thrown for unsupported locations. Setting '
            'test_throw=False...')
        test_throw = False

boutcore.init('-q -q -q -q')

results = []

if args.operator is not None:
    # just test the named operator
    try:
        arguments = single_argument_operators[args.operator]
        results += cycle_staggering(*arguments)
    except KeyError:
        try:
            arguments = two_argument_operators[args.operator]
            results += cycle_staggering2(*arguments)
        except KeyError:
            try:
                arguments = single_argument_derivative_operators[args.operator]
                results += cycle_staggering(*arguments)
            except KeyError:
                try:
                    arguments = two_argument_derivative_operators[args.operator]
                    results += cycle_staggering2(*arguments)
                except KeyError:
                    raise KeyError("Unrecognized operator: "+args.operator)
else:
    # Test all operators

    # single-argument operators
    for arguments in single_argument_operators.values():
        results += cycle_staggering(*arguments)

    # two-argument operators
    for arguments in two_argument_operators.values():
        results += cycle_staggering2(*arguments) # no staggering allowed

    if test_deriv_ops:
        # test derivative operators
        for arguments in single_argument_derivative_operators.values():
            results += cycle_staggering(*arguments)
        for arguments in two_argument_derivative_operators.values():
            results += cycle_staggering2(*arguments)

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
