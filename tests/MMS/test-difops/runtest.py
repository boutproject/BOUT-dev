#!/usr/bin/env python3

import boutcore
from boutdata.mms_alternate import *
import numpy
import sympy
import re
from copy import copy

class DifopsMMS:
    """
    Common methods/data for MMS testing of differential operators
    """
    def __init__(self, metric, fullTest=True, plotError=False):
        """
        Input metric is a Metric object as defined in boutdata.mms_alternate
        """

        self.testThrow = fullTest
        self.plotError = plotError

        self.metric = metric
        self.fullTest = fullTest
        self.meshDict = {}
        self.inputDict = {}
        self.inputDict2 = {}
        self.dimStaggerDict = {}
        self.dimStaggerExpectThrowDict = {}
        self.mxg = 2
        self.myg = 2

        self.min_exponent = 6
        self.max_exponent = 7
        self.ngrids = numpy.logspace(self.min_exponent, self.max_exponent, num=self.max_exponent-self.min_exponent+1, base=2).astype(int)

        # functions to use as input to operators
        self.testfunc = 'cos(2*pi*x+y+z)'
        self.testfunc2 = 'sin(4*pi*x+2*y+2*z)+cos(2*pi*x-z)'

        self.analytic_input = sympy.sympify(self.testfunc)
        self.analytic_input2 = sympy.sympify(self.testfunc2)

        self.locations = ['CENTRE', 'XLOW', 'YLOW', 'ZLOW']
        #self.ftypes = ['3D', '2D']
        self.ftypes = ['3D']

        # list of test results
        self.results = []

        # search pattern for splitKey
        self.keyRE = re.compile('(nx)([0-9]+)(ny)([0-9]+)(nz)([0-9]+)(.*)')

        if self.testThrow:
            if boutcore.bout_CHECK < 1:
                print('Warning: CHECK='+str(boutcore.bout_CHECK)+' so exceptions will '
                    'not be thrown for unsupported locations. Setting '
                    'testThrow=False...')
                self.testThrow = False
            else:
                # add small mesh and inputs to test for expected exceptions
                self.meshDict['expectThrow'] = self.makeMesh(64, 64, 3)
                for location in self.locations:
                    for ftype in self.ftypes:
                        self.inputDict['expectThrow'+location+ftype] = self.makeField(self.testfunc, 'expectThrow', location, ftype)

    def getKeyBase(self, n, dimensions, loc):
        """
        Given a grid size and strings defining the dimensions and staggering needed, create a key specifing nx, ny, nz
        """

        if 'x' in dimensions or loc=='XLOW':
            nx = n
        else:
            nx = 1
        if 'y' in dimensions or loc=='YLOW':
            ny = n
        else:
            ny = 1
        if 'z' in dimensions or loc=='ZLOW':
            nz = n
        else:
            nz = 1

        return 'nx'+str(nx)+'ny'+str(ny)+'nz'+str(nz)

    def splitKey(self, key):
        """
        Convert a string into {nx,ny,nz,stagger}
        """

        #if key[-1] in 'sxyz':
        #    stagger = key[-1]
        #    key = key[:-1]
        #else:
        #    stagger = False

        #if key[:2] == 'nx':
        #    key = key[2:]
        #else:
        #    throw ValueError("incorrect key format")

        #nx, key = key.split('ny')
        #nx = int(nx)

        #ny, key = key.split('nz')
        #ny = int(ny)

        #try:
        #    nz = int(key)
        #except:
        #    throw ValueError("Mesh size key parsing error: could not convert "+key+" to int for nz")

        result = self.keyRE.match(key)
        nx = int(result.group(2))
        ny = int(result.group(4))
        nz = int(result.group(6))
        stagger = result.group(7)

        return nx, ny, nz, stagger

    def getMesh(self, key):
        try:
            return self.meshDict[key]
        except KeyError:
            mesh = self.makeMesh(*self.splitKey(key))
            self.meshDict[key] = mesh
            return mesh

    def makeMesh(self, nx, ny, nz, stagger=''):
        if nx > 1:
            mxg = self.mxg
        else:
            mxg = 0
        if ny > 1:
            myg = self.myg
        else:
            myg = 0

        # set options
        boutcore.setOption('mxg', str(mxg), force=True)
        boutcore.setOption('myg', str(myg), force=True)
        # set up mesh input
        boutcore.setOption('testmesh:nx', exprToStr(nx+2*mxg), force=True)
        boutcore.setOption('testmesh:ny', exprToStr(ny), force=True)
        boutcore.setOption('testmesh:nz', exprToStr(nz), force=True)
        boutcore.setOption('testmesh:dx', exprToStr(self.metric.psiwidth*self.metric.scalex/nx), force=True)
        boutcore.setOption('testmesh:dy', exprToStr(2.*pi*self.metric.scaley/ny), force=True)
        boutcore.setOption('testmesh:dz', exprToStr(2.*pi/nz), force=True)
        boutcore.setOption('testmesh:g11', exprToStr(self.metric.g11), force=True)
        boutcore.setOption('testmesh:g22', exprToStr(self.metric.g22), force=True)
        boutcore.setOption('testmesh:g33', exprToStr(self.metric.g33), force=True)
        boutcore.setOption('testmesh:g12', exprToStr(self.metric.g12), force=True)
        boutcore.setOption('testmesh:g13', exprToStr(self.metric.g13), force=True)
        boutcore.setOption('testmesh:g23', exprToStr(self.metric.g23), force=True)
        boutcore.setOption('testmesh:g_11', exprToStr(self.metric.g_11), force=True)
        boutcore.setOption('testmesh:g_22', exprToStr(self.metric.g_22), force=True)
        boutcore.setOption('testmesh:g_33', exprToStr(self.metric.g_33), force=True)
        boutcore.setOption('testmesh:g_12', exprToStr(self.metric.g_12), force=True)
        boutcore.setOption('testmesh:g_13', exprToStr(self.metric.g_13), force=True)
        boutcore.setOption('testmesh:g_23', exprToStr(self.metric.g_23), force=True)
        boutcore.setOption('testmesh:staggergrids', str('true'), force=True)

        # create new Mesh object
        return boutcore.Mesh(section='testmesh')

    def getInput(self, keyBase, loc, ftype):
        try:
            return self.inputDict[keyBase+loc+ftype]
        except KeyError:
            inputField = self.makeField(self.testfunc, keyBase, loc, ftype)
            self.inputDict[keyBase+loc+ftype] = inputField
            return inputField

    def getInput2(self, key, ftype):
        try:
            return self.inputDict2[key+loc+ftype]
        except KeyError:
            inputField = self.makeField(self.testfunc2, keyBase, loc, ftype)
            self.inputDict2[key+loc+ftype] = inputField
            return inputField

    def makeField(self, func, keyBase, loc, ftype):
        # ensure we don't change the input 'func' string
        func = copy(func)
        mesh = self.getMesh(keyBase)

        if ftype == '2D':
            # cannot have z-dependence
            func = func.replace('z', '0')
            return boutcore.create2D(func, mesh, outloc=loc)
        elif ftype == '3D':
            return boutcore.create3D(func, mesh, outloc=loc)
        else:
            raise ValueError('Unexpected ftype argument '+str(ftype))

    def getDimStagger(self, base_dimensions, stagger_directions):
        try:
            return self.dimStaggerDict[base_dimensions][stagger_directions]
        except KeyError:
            # all derivatives at same inloc/outloc should work
            dimensions_staggers = [(base_dimensions, ('CENTRE', 'CENTRE')),
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
            try:
                self.dimStaggerDict[base_dimensions][stagger_directions] = dimensions_staggers
            except KeyError:
                self.dimStaggerDict[base_dimensions] = {}
                self.dimStaggerDict[base_dimensions][stagger_directions] = dimensions_staggers

            return dimensions_staggers

    def getDimStaggerExpectThrow(self, base_dimensions, stagger_directions):
        try:
            return self.dimStaggerExpectThrowDict[base_dimensions][stagger_directions]
        except KeyError:
            # first make a list of all permutations
            fail_staggers = [(x,y) for x in self.locations for y in self.locations]
            dimensions_staggers = self.getDimStagger(base_dimensions, stagger_directions)
            for dimensions, stagger in dimensions_staggers:
                if stagger is not None:
                    index = fail_staggers.index(stagger)
                    del fail_staggers[index]
            fail_dimensions_staggers = [('xyz', s) for s in fail_staggers]
            try:
                self.dimStaggerExpectThrowDict[base_dimensions][stagger_directions] = fail_dimensions_staggers
            except KeyError:
                self.dimStaggerExpectThrowDict[base_dimensions] = {}
                self.dimStaggerExpectThrowDict[base_dimensions][stagger_directions] = fail_dimensions_staggers
            return fail_dimensions_staggers

    def testOperatorAtLocation(self, dimensions, boutcore_operator, symbolic_operator, order, ftype, method, stagger):
        error_list = []
        if full_test:
            print('testing',boutcore_operator, ftype, stagger)
        inloc = stagger[0]
        outloc = stagger[1]

        # stag records if either of the locations is staggered
        if inloc == 'CENTRE':
            stag = outloc
        else:
            stag = inloc

        print('testing',boutcore_operator, ftype, stagger)
        error_list = []
        for n in self.ngrids:
            if full_test:
                print('n =',n)
            keyBase = self.getKeyBase(n, dimensions, stag)
            mesh = self.getMesh(keyBase)

            # calculate result of differential operator using BOUT++ implementation
            bout_input = self.getInput(keyBase, inloc, ftype)
            if method is None:
                bout_result = boutcore_operator(bout_input, outloc=outloc)
            else:
                bout_result = boutcore_operator(bout_input, outloc=outloc, method=method)

            # calculate result of differential operator symbolically, then convert to boutcore.Field3D/Field2D
            analytic_func = symbolic_operator(self.analytic_input)
            if ftype == '2D':
                analytic_result = boutcore.create2D(exprToStr(analytic_func), mesh, outloc=outloc)
            elif ftype == '3D':
                analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh, outloc=outloc)
            else:
                raise ValueError('Unexpected ftype argument '+str(ftype))

            # calculate max error
            error = bout_result - analytic_result # as Field3D/Field2D
            error = error.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1] # numpy array, without guard cells
            error_list.append(numpy.max(numpy.abs(error))) # max error

        logerrors = numpy.log(error_list[-2]/error_list[-1])
        logspacing = numpy.log(self.ngrids[-1]/self.ngrids[-2])
        convergence = logerrors/logspacing

        if order-.1 < convergence < order+.2:
            return 'pass'
        else:
            if self.plotError:
                from matplotlib import pyplot
                pyplot.loglog(1./ngrids, error_list)
                pyplot.show()
                from boututils.showdata import showdata
                showdata(error)
            return str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' '+str(ftype)+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.'

    def expectThrowAtLocation(self, boutcore_operator, ftype, method, stagger):
        inloc = stagger[0]
        outloc = stagger[1]

        print('testing',boutcore_operator, ftype, stagger)

        try:
            boutcore_operator(self.getInput('expectThrow', inloc, ftype), outloc=outloc)
        except RuntimeError:
            return 'pass'
        else:
            return 'Expected '+str(boutcore_operator)+' to throw for '+stagger[0]+'->'+stagger[1]+' '+' '+str(ftype)+' '+str(method)+' but it did not.'

    def testOperator(self, stagger_directions, base_dimensions, boutcore_operator, symbolic_operator, order, ftype, method=None):
        for dimensions,stagger in self.getDimStagger(base_dimensions, stagger_directions):
            self.results.append(self.testOperatorAtLocation(dimensions, boutcore_operator, symbolic_operator, order, ftype, method, stagger))

        if self.testThrow:
            for dimensions,stagger in self.getDimStaggerExpectThrow(base_dimensions, stagger_directions):
                self.results.append(self.expectThrowAtLocation(boutcore_operator, ftype, method, stagger))

    def checkResults(self):
        fail = False
        for result in self.results:
            if result is not 'pass':
                print(result)
                fail = True
        return not fail


if __name__ == "__main__":
    from sys import exit

    # get command line arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--short', action='store_true', default=False)
    parser.add_argument('--plot', action='store_true', default=False)
    args = parser.parse_args()
    full_test = not args.short

    # geometry for simple circular tokamak
    tokamak = SimpleTokamak()

    # rescale x and y coordinates so dx and dy are not constants
    tokamak.set_scalex(1 + .1*sin(2*pi*metric.x+metric.y))
    tokamak.set_scaley(1 + .1*sin(2*pi*metric.x-metric.y))
    # re-calculate metric terms
    tokamak.metric()

    boutcore.init('-q -q -q -q')

    driver = DifopsMMS(metric, full_test, args.plot)

    driver.testOperator('y', 'y', boutcore.Grad_par, Grad_par, 2, '3D')

    if driver.checkResults():
        print('pass')
        exit(0)
    else:
        print('fail')
        exit(1)
