#!/usr/bin/env python3

import boutcore
from boutdata.mms_alternate import *
import numpy
import sympy
import re
from copy import copy
from collections import OrderedDict

class DifopsMMS:
    """
    Common methods/data for MMS testing of differential operators
    """
    def __init__(self, metric, fullTest=True, test3D=True, plotError=False):
        """
        Input metric is a Metric object as defined in boutdata.mms_alternate
        """

        self.metric = metric
        self.fullTest = fullTest
        self.testThrow = fullTest
        self.test3D = test3D
        self.plotError = plotError
        self.meshDict = {}
        self.inputDict = {}
        self.dimStaggerDict = {}
        self.dimStaggerExpectThrowDict = {}
        self.inputDict2 = {}
        self.dimStaggerDict2 = {}
        self.dimStaggerExpectThrowDict2 = {}
        self.mxg = 2
        self.myg = 2

        if self.fullTest:
            # Use high enough resolution for all tests to reach expected order
            self.min_exponent = 8
            self.max_exponent = 9
        else:
            # run at low resolution, and allow faster than expected convergence
            # for some tests where this is expected behaviour
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
                self.meshDict['expectThrow'] = self.makeMesh(2, 12, 3)
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

    def getInput2(self, keyBase, loc, ftype):
        try:
            return self.inputDict2[keyBase+loc+ftype]
        except KeyError:
            inputField = self.makeField(self.testfunc2, keyBase, loc, ftype)
            self.inputDict2[keyBase+loc+ftype] = inputField
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

            # now remove the ones that should NOT throw
            dimensions_staggers = self.getDimStagger(base_dimensions, stagger_directions)
            for dimensions, stagger in dimensions_staggers:
                index = fail_staggers.index(stagger)
                del fail_staggers[index]

            fail_dimensions_staggers = [('xyz', s) for s in fail_staggers]
            try:
                self.dimStaggerExpectThrowDict[base_dimensions][stagger_directions] = fail_dimensions_staggers
            except KeyError:
                self.dimStaggerExpectThrowDict[base_dimensions] = {}
                self.dimStaggerExpectThrowDict[base_dimensions][stagger_directions] = fail_dimensions_staggers
            return fail_dimensions_staggers

    def testOperatorAtLocation(self, dimensions, boutcore_operator, symbolic_operator, order, ftype, fudged_max_order, method, stagger):
        error_list = []
        print('testing',boutcore_operator, ftype, stagger)
        inloc = stagger[0]
        outloc = stagger[1]

        # stag records if either of the locations is staggered
        if inloc == 'CENTRE':
            stag = outloc
        else:
            stag = inloc

        if ftype == '2D':
            # cannot have z-dependence
            analytic_input = self.analytic_input.replace('z', '0')
        else:
            analytic_input = self.analytic_input

        error_list = []
        errors = []
        boutfields = []
        sympyfields = []
        for n in self.ngrids:
            if self.fullTest:
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
            analytic_func = symbolic_operator(analytic_input)
            if ftype == '2D':
                analytic_result = boutcore.create2D(exprToStr(analytic_func), mesh, outloc=outloc)
            elif ftype == '3D':
                analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh, outloc=outloc)
            else:
                raise ValueError('Unexpected ftype argument '+str(ftype))

            # calculate max error
            error = bout_result - analytic_result # as Field3D/Field2D
            save_inds = numpy.index_exp[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1]
            errors.append(error.get()[save_inds])
            boutfields.append(bout_result.get()[save_inds])
            sympyfields.append(analytic_result.get()[save_inds])
            error = error.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1] # numpy array, without guard cells
            error_list.append(numpy.max(numpy.abs(error))) # max error

        logerrors = numpy.log(error_list[-2]/error_list[-1])
        logspacing = numpy.log(self.ngrids[-1]/self.ngrids[-2])
        convergence = logerrors/logspacing

        if self.fullTest or fudged_max_order is None:
            # full, strict test for convergence order
            max_order = order+.2
        else:
            # allow some tests to expect a higher order of convergence for
            # quick tests at lower resolution
            max_order = fudged_max_order

        if order-.1 < convergence < max_order:
            return 'pass'
        else:
            error_string = str(boutcore_operator)+' is not working for '+inloc+'->'+outloc+' '+str(ftype)+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.'
            if self.plotError:
                print(error_string)
                from matplotlib import pyplot
                pyplot.loglog(1./self.ngrids, error_list, label="sim")
                pyplot.loglog(1./self.ngrids, self.ngrids[-1]**order/self.ngrids**order*error_list[-1], 'k--', label="expected order")
                pyplot.legend()
                pyplot.show()
                from boututils.showdata import showdata
                plot_error = copy(error)
                plot_error = numpy.squeeze(plot_error)
                pb = numpy.squeeze(bout_result.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1])
                ps = numpy.squeeze(analytic_result.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1])
                if len(plot_error.shape) == 1:
                    plot_error = plot_error[numpy.newaxis, :]
                    pb = pb[numpy.newaxis, :]
                    ps = ps[numpy.newaxis, :]
                #showdata(plot_error)
                showdata([plot_error,pb,ps],titles=['error','bout','sympy'])
                #showdata([plot_error, numpy.squeeze(mesh.getCoordinates().g_22.get())[None,:], numpy.squeeze(mesh.getCoordinates('YLOW').g_22.get())[None,:]])
                #for f in [(metric.g_22, mesh.getCoordinates().g_22)]:
                #for fa,fb,name in [(metric.g11, mesh.getCoordinates().g11, 'g11'),
                #                   (metric.g22, mesh.getCoordinates().g22, 'g22'),
                #                   (metric.g33, mesh.getCoordinates().g33, 'g33'),
                #                   (metric.g12, mesh.getCoordinates().g12, 'g12'),
                #                   (metric.g13, mesh.getCoordinates().g13, 'g13'),
                #                   (metric.g23, mesh.getCoordinates().g23, 'g23'),
                #                   (metric.G1, mesh.getCoordinates().G1, 'G1'),
                #                   (metric.G2, mesh.getCoordinates().G2, 'G2'),
                #                   (metric.G3, mesh.getCoordinates().G3, 'G3')
                #                  ]:
                #    fa=boutcore.create2D(exprToStr(fa), mesh, outloc=outloc)
                #    fa=numpy.squeeze(fa.get()[mesh.xstart:mesh.xend+1,mesh.ystart:mesh.yend+1])[None,:]
                #    fb=numpy.squeeze(fb.get()[mesh.xstart:mesh.xend+1,mesh.ystart:mesh.yend+1])[None,:]
                #    print(fa.shape, fb.shape)
                #    showdata([fb, fa, fb-fa], titles=[name+' bout',name+' sympy',name+' diff'])
                pyplot.figure()
                for e,b,s in zip(errors,boutfields,sympyfields):
                    inds = numpy.index_exp[:,2]
                    e = e[inds]
                    b = b[inds]
                    s = s[inds]
                    x = numpy.linspace(0.,1.,e.shape[0])
                    pyplot.subplot(131)
                    pyplot.semilogy(x,numpy.abs(e), label=x.shape[0])
                    pyplot.subplot(132)
                    pyplot.plot(x, b, label=x.shape[0])
                    pyplot.subplot(133)
                    pyplot.plot(x, s, label=x.shape[0])
                pyplot.title(str(boutcore_operator))
                pyplot.subplot(131)
                pyplot.title('error')
                pyplot.legend()
                pyplot.subplot(132)
                pyplot.title('bout')
                pyplot.legend()
                pyplot.subplot(133)
                pyplot.title('sympy')
                pyplot.legend()
                pyplot.show()
            return error_string

    def expectThrowAtLocation(self, boutcore_operator, ftype, method, stagger):
        inloc = stagger[0]
        outloc = stagger[1]

        print('testing',boutcore_operator, ftype, stagger)

        try:
            if method is None:
                boutcore_operator(self.getInput('expectThrow', inloc, ftype), outloc=outloc)
            else:
                boutcore_operator(self.getInput('expectThrow', inloc, ftype), outloc=outloc, method=method)
        except RuntimeError:
            return 'pass'
        else:
            return 'Expected '+str(boutcore_operator)+' to throw for '+stagger[0]+'->'+stagger[1]+' '+' '+str(ftype)+' '+str(method)+' but it did not.'

    def testOperator(self, stagger_directions, base_dimensions, boutcore_operator, symbolic_operator, order, ftype, fudged_max_order=None, method=None):
        for dimensions,stagger in self.getDimStagger(base_dimensions, stagger_directions):
            if (not self.test3D) and 'x' in dimensions and 'y' in dimensions and 'z' in dimensions:
                continue
            self.results.append(self.testOperatorAtLocation(dimensions, boutcore_operator, symbolic_operator, order, ftype, fudged_max_order, method, stagger))

        if self.testThrow:
            for dimensions,stagger in self.getDimStaggerExpectThrow(base_dimensions, stagger_directions):
                self.results.append(self.expectThrowAtLocation(boutcore_operator, ftype, method, stagger))

    def getDimStagger2(self, base_dimensions, stagger_directions):
        try:
            return self.dimStaggerDict2[base_dimensions][stagger_directions]
        except KeyError:
            # all derivatives at same inloc/outloc should work
            dimensions_staggers = [(base_dimensions, ('CENTRE', 'CENTRE', 'CENTRE')),
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
            try:
                self.dimStaggerDict2[base_dimensions][stagger_directions] = dimensions_staggers
            except KeyError:
                self.dimStaggerDict2[base_dimensions] = {}
                self.dimStaggerDict2[base_dimensions][stagger_directions] = dimensions_staggers

            return dimensions_staggers

    def getDimStaggerExpectThrow2(self, base_dimensions, stagger_directions):
        try:
            return self.dimStaggerExpectThrowDict2[base_dimensions][stagger_directions]
        except KeyError:
            # first make a list of all permutations
            fail_staggers = [(x,y,z) for x in self.locations for y in self.locations for z in self.locations]

            # now remove the ones that should NOT throw
            dimensions_staggers = self.getDimStagger2(base_dimensions, stagger_directions)
            for dimensions, stagger in dimensions_staggers:
                index = fail_staggers.index(stagger)
                del fail_staggers[index]

            fail_dimensions_staggers = [('xyz', s) for s in fail_staggers]
            try:
                self.dimStaggerExpectThrowDict2[base_dimensions][stagger_directions] = fail_dimensions_staggers
            except KeyError:
                self.dimStaggerExpectThrowDict2[base_dimensions] = {}
                self.dimStaggerExpectThrowDict2[base_dimensions][stagger_directions] = fail_dimensions_staggers
            return fail_dimensions_staggers

    def testOperatorAtLocation2(self, dimensions, boutcore_operator, symbolic_operator, order, ftype, fudged_max_order, method, stagger):
        error_list = []
        print('testing',boutcore_operator, ftype, stagger)
        inloc1 = stagger[0]
        inloc2 = stagger[1]
        outloc = stagger[2]

        # stag records if any of the locations is staggered
        stag = 'CENTRE'
        if inloc1 != 'CENTRE':
            stag = inloc1
        elif inloc2 != 'CENTRE':
            stag = inloc2
        elif outloc != 'CENTRE':
            stag = outloc

        if ftype == '2D':
            # cannot have z-dependence
            analytic_input1 = self.analytic_input.replace('z', '0')
            analytic_input2 = self.analytic_input2.replace('z', '0')
        else:
            analytic_input1 = self.analytic_input
            analytic_input2 = self.analytic_input2

        error_list = []
        errors = []
        boutfields = []
        sympyfields = []
        for n in self.ngrids:
            if self.fullTest:
                print('n =',n)
            keyBase = self.getKeyBase(n, dimensions, stag)
            mesh = self.getMesh(keyBase)

            # calculate result of differential operator using BOUT++ implementation
            bout_input1 = self.getInput(keyBase, inloc1, ftype)
            bout_input2 = self.getInput2(keyBase, inloc2, ftype)
            if method is None:
                bout_result = boutcore_operator(bout_input1, bout_input2, outloc=outloc)
            else:
                bout_result = boutcore_operator(bout_input1, bout_input2, outloc=outloc, method=method)

            # calculate result of differential operator symbolically, then convert to boutcore.Field3D/Field2D
            analytic_func = symbolic_operator(analytic_input1, analytic_input2)
            if ftype == '2D':
                analytic_result = boutcore.create2D(exprToStr(analytic_func), mesh, outloc=outloc)
            elif ftype == '3D':
                analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh, outloc=outloc)
            else:
                raise ValueError('Unexpected ftype argument '+str(ftype))

            # calculate max error
            error = bout_result - analytic_result # as Field3D/Field2D
            save_inds = numpy.index_exp[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1]
            errors.append(error.get()[save_inds])
            boutfields.append(bout_result.get()[save_inds])
            sympyfields.append(analytic_result.get()[save_inds])
            error = error.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1] # numpy array, without guard cells
            error_list.append(numpy.max(numpy.abs(error))) # max error

        logerrors = numpy.log(error_list[-2]/error_list[-1])
        logspacing = numpy.log(self.ngrids[-1]/self.ngrids[-2])
        convergence = logerrors/logspacing

        if self.fullTest or fudged_max_order is None:
            # full, strict test for convergence order
            max_order = order+.2
        else:
            # allow some tests to expect a higher order of convergence for
            # quick tests at lower resolution
            max_order = fudged_max_order

        if order-.1 < convergence < max_order:
            return 'pass'
        else:
            error_string = str(boutcore_operator)+' is not working for {'+inloc1+','+inloc2+'}->'+outloc+' '+str(ftype)+' '+str(method)+'. Expected '+str(order)+', got '+str(convergence)+'.'
            if self.plotError:
                print(error_string)
                from matplotlib import pyplot
                pyplot.loglog(1./self.ngrids, error_list, label="sim")
                pyplot.loglog(1./self.ngrids, self.ngrids[-1]**order/self.ngrids**order*error_list[-1], 'k--', label="expected order")
                pyplot.legend()
                pyplot.show()
                from boututils.showdata import showdata
                plot_error = copy(error)
                plot_error = numpy.squeeze(plot_error)
                pb = numpy.squeeze(bout_result.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1])
                ps = numpy.squeeze(analytic_result.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1])
                if len(plot_error.shape) == 1:
                    plot_error = plot_error[numpy.newaxis, :]
                    pb = pb[numpy.newaxis, :]
                    ps = ps[numpy.newaxis, :]
                #showdata(plot_error)
                showdata([plot_error,pb,ps],titles=['error','bout','sympy'])
                pyplot.figure()
                for e,b,s in zip(errors,boutfields,sympyfields):
                    inds = numpy.index_exp[:,2]
                    e = e[inds]
                    b = b[inds]
                    s = s[inds]
                    x = numpy.linspace(0.,1.,e.shape[0])
                    pyplot.subplot(131)
                    pyplot.semilogy(x,numpy.abs(e), label=x.shape[0])
                    pyplot.subplot(132)
                    pyplot.plot(x, b, label=x.shape[0])
                    pyplot.subplot(133)
                    pyplot.plot(x, s, label=x.shape[0])
                pyplot.title(str(boutcore_operator))
                pyplot.subplot(131)
                pyplot.title('error')
                pyplot.legend()
                pyplot.subplot(132)
                pyplot.title('bout')
                pyplot.legend()
                pyplot.subplot(133)
                pyplot.title('sympy')
                pyplot.legend()
                pyplot.show()
            return error_string

    def expectThrowAtLocation2(self, boutcore_operator, ftype, method, stagger):
        inloc1 = stagger[0]
        inloc2 = stagger[1]
        outloc = stagger[2]

        print('testing',boutcore_operator, ftype, stagger)

        try:
            if method is None:
                boutcore_operator(self.getInput('expectThrow', inloc1, ftype), self.getInput2('expectThrow', inloc2, ftype), outloc=outloc)
            else:
                boutcore_operator(self.getInput('expectThrow', inloc1, ftype), self.getInput2('expectThrow', inloc2, ftype), outloc=outloc, method=method)
        except RuntimeError:
            return 'pass'
        else:
            return 'Expected '+str(boutcore_operator)+' to throw for {'+inloc1+','+inloc2+'}->'+outloc+' '+' '+str(ftype)+' '+str(method)+' but it did not.'

    def testOperator2(self, stagger_directions, base_dimensions, boutcore_operator, symbolic_operator, order, ftype, fudged_max_order=None, method=None):
        for dimensions,stagger in self.getDimStagger2(base_dimensions, stagger_directions):
            if (not self.test3D) and 'x' in dimensions and 'y' in dimensions and 'z' in dimensions:
                continue
            self.results.append(self.testOperatorAtLocation2(dimensions, boutcore_operator, symbolic_operator, order, ftype, fudged_max_order, method, stagger))

        if self.testThrow:
            for dimensions,stagger in self.getDimStaggerExpectThrow2(base_dimensions, stagger_directions):
                self.results.append(self.expectThrowAtLocation2(boutcore_operator, ftype, method, stagger))

    def checkResults(self):
        fail = False
        for result in self.results:
            if result is not 'pass':
                print(result)
                fail = True
        return not fail

    def testField(self, boutcore_field, symbolic_field, order, ftype, stagger):
        error_list = []
        print('testing',boutcore_field, ftype, stagger)

        error_list = []
        errors = []
        boutfields = []
        sympyfields = []
        for n in self.ngrids:
            if self.fullTest:
                print('n =',n)
            keyBase = self.getKeyBase(n, 'xy', stagger)
            mesh = self.getMesh(keyBase)

            # result using BOUT++ implementation
            if boutcore_field == 'g11':
                bout_result = mesh.getCoordinates(stagger).g11
            elif boutcore_field == 'g22':
                bout_result = mesh.getCoordinates(stagger).g22
            elif boutcore_field == 'g33':
                bout_result = mesh.getCoordinates(stagger).g33
            elif boutcore_field == 'g12':
                bout_result = mesh.getCoordinates(stagger).g12
            elif boutcore_field == 'g13':
                bout_result = mesh.getCoordinates(stagger).g13
            elif boutcore_field == 'g23':
                bout_result = mesh.getCoordinates(stagger).g23
            elif boutcore_field == 'dx':
                bout_result = mesh.getCoordinates(stagger).dx
            elif boutcore_field == 'dy':
                bout_result = mesh.getCoordinates(stagger).dy
            elif boutcore_field == 'J':
                bout_result = mesh.getCoordinates(stagger).J
            elif boutcore_field == 'G1':
                bout_result = mesh.getCoordinates(stagger).G1
            elif boutcore_field == 'G2':
                bout_result = mesh.getCoordinates(stagger).G2
            elif boutcore_field == 'G3':
                bout_result = mesh.getCoordinates(stagger).G3
            else:
                raise ValueError("unsupported field")

            # result using symbolic input, then convert to boutcore.Field3D/Field2D
            if ftype == '2D':
                analytic_result = boutcore.create2D(exprToStr(symbolic_field), mesh, outloc=stagger)
            elif ftype == '3D':
                analytic_result = boutcore.create3D(exprToStr(analytic_func), mesh, outloc=stagger)
            else:
                raise ValueError('Unexpected ftype argument '+str(ftype))

            # calculate max error
            error = bout_result - analytic_result # as Field3D/Field2D
            save_inds = numpy.index_exp[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1]
            errors.append(error.get()[save_inds])
            boutfields.append(bout_result.get()[save_inds])
            sympyfields.append(analytic_result.get()[save_inds])
            error = error.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1] # numpy array, without guard cells
            error_list.append(numpy.max(numpy.abs(error))) # max error
            #errors.append(error)

        logerrors = numpy.log(error_list[-2]/error_list[-1])
        logspacing = numpy.log(self.ngrids[-1]/self.ngrids[-2])
        convergence = logerrors/logspacing

        if order-.1 < convergence < order+.2:
            return 'pass'
        else:
            error_string = str(boutcore_field)+' is not working for '+stagger+' '+str(ftype)+'. Expected '+str(order)+', got '+str(convergence)+'.'
            if self.plotError:
                print(error_string)
                from matplotlib import pyplot
                pyplot.loglog(1./self.ngrids, error_list, label="sim")
                pyplot.loglog(1./self.ngrids, self.ngrids[-1]**order/self.ngrids**order*error_list[-1], 'k--', label="expected order")
                pyplot.legend()
                pyplot.show()
                from boututils.showdata import showdata
                #plot_error = copy(error)
                #plot_error = numpy.squeeze(plot_error)
                #pb = numpy.squeeze(bout_result.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1])
                #ps = numpy.squeeze(analytic_result.get()[mesh.xstart:mesh.xend+1, mesh.ystart:mesh.yend+1])
                #if len(plot_error.shape) == 1:
                #    plot_error = plot_error[numpy.newaxis, :]
                #    pb = pb[numpy.newaxis, :]
                #    ps = ps[numpy.newaxis, :]
                #else:
                #    plot_error = plot_error[:3]
                #    pb = pb[:3]
                #    ps = ps[:3]
                ##showdata(plot_error)
                #showdata([plot_error,pb,ps],titles=['error','bout','sympy'])
                pyplot.figure()
                for e,b,s in zip(errors,boutfields,sympyfields):
                    inds = numpy.index_exp[:,2]
                    e = e[inds]
                    b = b[inds]
                    s = s[inds]
                    x = numpy.linspace(0.,1.,e.shape[0])
                    pyplot.subplot(131)
                    pyplot.semilogy(x,numpy.abs(e), label=x.shape[0])
                    pyplot.subplot(132)
                    pyplot.plot(x, b, label=x.shape[0])
                    pyplot.subplot(133)
                    pyplot.plot(x, s, label=x.shape[0])
                pyplot.title(boutcore_field)
                pyplot.subplot(131)
                pyplot.title('error')
                pyplot.legend()
                pyplot.subplot(132)
                pyplot.title('bout')
                pyplot.legend()
                pyplot.subplot(133)
                pyplot.title('sympy')
                pyplot.legend()
                pyplot.show()
            return error_string

if __name__ == "__main__":
    from sys import exit

    # get command line arguments
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--short', action='store_true', default=False)
    parser.add_argument('--test3D', action='store_true', default=False)
    parser.add_argument('--plot', action='store_true', default=False)
    parser.add_argument('--operator', default=None)
    parser.add_argument('-v', action='store_true', default=False)
    args = parser.parse_args()
    full_test = not args.short

    # geometry for simple circular tokamak
    geometry = PeriodicGeometry()

    # rescale x and y coordinates so dx and dy are not constants
    geometry.set_scalex(1 + .1*sin(2*pi*metric.x+metric.y))
    geometry.set_scaley(1 + .1*sin(2*pi*metric.x-metric.y))
    # re-calculate metric terms
    geometry.metric()

    if args.v:
        # verbose output
        boutcore.init('-v -v')
    else:
        boutcore.init('-q -q -q -q')

    # create the testing class
    driver = DifopsMMS(metric, full_test, args.test3D, args.plot)

    # Store the inputs in a dict so we can look them up or iterate through them.
    # Note: Several operators on Field2D give exactly zero, so they have no
    # discretization error and we cannot test them with MMS. They are therefore
    # excluded from this script.
    operator_inputs = OrderedDict([
            ( 'Grad_par', ('y', 'y', boutcore.Grad_par, Grad_par, 2, '3D') ),
            ( 'Div_par', ('y', 'y', boutcore.Div_par, Div_par, 2, '3D', 3.6) ),
            ( 'Grad2_par2', ('y', 'y', boutcore.Grad2_par2, Grad2_par2, 2, '3D') ),
            ( 'Laplace', ('', 'xyz', boutcore.Laplace, Laplace, 2, '3D') ),
            ( 'Laplace_par', ('y', 'y', boutcore.Laplace_par, Laplace_par, 2, '3D', 2.5) ),
            ( 'Laplace_perp', ('', 'xyz', boutcore.Laplace_perp, Laplace_perp, 2, '3D') ),
            ( 'DDX', ('x', 'x', boutcore.DDX, DDX, 2, '3D') ),
            ( 'DDY', ('y', 'y', boutcore.DDY, DDY, 2, '3D') ),
            ( 'DDZ', ('z', 'z', boutcore.DDZ, DDZ, 2, '3D') ),
            ( 'D2DX2', ('x', 'x', boutcore.D2DX2, D2DX2, 2, '3D') ),
            ( 'D2DY2', ('y', 'y', boutcore.D2DY2, D2DY2, 2, '3D') ),
            ( 'D2DZ2', ('z', 'z', boutcore.D2DZ2, D2DZ2, 2, '3D') ),
            ( 'D4DX4', ('', 'x', boutcore.D4DX4, D4DX4, 2, '3D') ),
            ( 'D4DY4', ('', 'y', boutcore.D4DY4, D4DY4, 2, '3D') ),
            ( 'D4DZ4', ('', 'z', boutcore.D4DZ4, D4DZ4, 2, '3D') ),
            ( 'D2DXDY', ('y', 'xy', boutcore.D2DXDY, D2DXDY, 2, '3D') ),
            ( 'D2DXDZ', ('', 'xz', boutcore.D2DXDZ, D2DXDZ, 2, '3D') ),
            ( 'D2DYDZ', ('', 'yz', boutcore.D2DYDZ, D2DYDZ, 2, '3D') ),
            ( 'Grad_par_2D', ('y', 'y', boutcore.Grad_par, Grad_par, 2, '2D') ),
            ( 'Div_par_2D', ('y', 'y', boutcore.Div_par, Div_par, 2, '2D', 3.6) ),
            ( 'Grad2_par2_2D', ('y', 'y', boutcore.Grad2_par2, Grad2_par2, 2, '2D') ),
            ( 'Laplace_2D', ('', 'xy', boutcore.Laplace, Laplace, 2, '2D') ),
            ( 'Laplace_par_2D', ('y', 'y', boutcore.Laplace_par, Laplace_par, 2, '2D', 2.5) ),
            ( 'Laplace_perp_2D', ('', 'xy', boutcore.Laplace_perp, Laplace_perp, 2, '2D') ),
            ( 'DDX_2D', ('x', 'x', boutcore.DDX, DDX, 2, '2D') ),
            ( 'DDY_2D', ('y', 'y', boutcore.DDY, DDY, 2, '2D') ),
            ( 'D2DX2_2D', ('x', 'x', boutcore.D2DX2, D2DX2, 2, '2D') ),
            ( 'D2DY2_2D', ('y', 'y', boutcore.D2DY2, D2DY2, 2, '2D') ),
            ( 'D4DX4_2D', ('', 'x', boutcore.D4DX4, D4DX4, 2, '2D') ),
            ( 'D4DY4_2D', ('', 'y', boutcore.D4DY4, D4DY4, 2, '2D') ),
            ( 'D2DXDY_2D', ('y', 'xy', boutcore.D2DXDY, D2DXDY, 2, '2D') ),
            ])

    operator_inputs2 = OrderedDict([
            ( 'Vpar_Grad_par', ('y', 'y', boutcore.Vpar_Grad_par, Vpar_Grad_par, 2, '3D') ),
            ( 'Div_par_K_Grad_par', ('yy', 'y', boutcore.Div_par_K_Grad_par, Div_par_K_Grad_par, 2, '3D') ),
            ( 'Div_par_flux', ('y', 'y', boutcore.Div_par_flux, lambda v,f: Div_par(v*f), 2, '3D') ), # uses first-order upwind method
            ( 'Div_par_flux_C2', ('', 'y', boutcore.Div_par_flux, lambda v,f: Div_par(v*f), 2, '3D', None, 'C2') ), # centred differencing
            ( 'bracket', ('', 'xz', boutcore.bracket, lambda f,g: bracket(f, g, include_yderivs=False), 2, '3D', None, 'BRACKET_ARAKAWA') ),
            ( 'bracket_STD', ('', 'xyz', boutcore.bracket, bracket, 2, '3D', 2.7, 'BRACKET_STD') ),
            ( 'VDDX', ('x', 'x', boutcore.VDDX, lambda v,f: v*DDX(f), 2, '3D') ),
            ( 'VDDY', ('y', 'y', boutcore.VDDY, lambda v,f: v*DDY(f), 2, '3D') ),
            ( 'VDDZ', ('z', 'z', boutcore.VDDZ, lambda v,f: v*DDZ(f), 2, '3D') ),
            ( 'FDDX', ('x', 'x', boutcore.FDDX, lambda v,f: DDX(v*f), 2, '3D') ),
            ( 'FDDY', ('y', 'y', boutcore.FDDY, lambda v,f: DDY(v*f), 2, '3D') ),
            ( 'FDDZ', ('z', 'z', boutcore.FDDZ, lambda v,f: DDZ(v*f), 2, '3D') ),
            ( 'Vpar_Grad_par_2D', ('y', 'y', boutcore.Vpar_Grad_par, Vpar_Grad_par, 2, '2D') ),
            ( 'bracket_STD_2D', ('', 'xyz', boutcore.bracket, bracket, 2, '2D', 2.8, 'BRACKET_STD') ),
            ( 'VDDX_2D', ('x', 'x', boutcore.VDDX, lambda v,f: v*DDX(f), 2, '2D') ),
            ( 'VDDY_2D', ('y', 'y', boutcore.VDDY, lambda v,f: v*DDY(f), 2, '2D') ),
            ( 'FDDX_2D', ('x', 'x', boutcore.FDDX, lambda v,f: DDX(v*f), 2, '2D') ),
            ( 'FDDY_2D', ('y', 'y', boutcore.FDDY, lambda v,f: DDY(v*f), 2, '2D') ),
            ])

    if args.operator in ['CENTRE', 'XLOW', 'YLOW', 'ZLOW']:
        stag = args.operator
        if stag == 'ZLOW':
            raise ValueError('No interpolation of Field2D for ZLOW, so does not make sense to check convergence')
        driver.testField('g11', metric.g11, 3, '2D', stag)
        # constant, so BOUT++ field is exact, so don't test #driver.testField('g22', metric.g22, 3, '2D', stag)
        driver.testField('g33', metric.g33, 3, '2D', stag)
        # zero, so don't test #driver.testField('g12', metric.g12, 3, '2D', stag)
        driver.testField('g13', metric.g13, 3, '2D', stag)
        if stag != 'XLOW':
            # quadratic in x, does not show convergence for XLOW as only
            # floating-point errors in XLOW interpolation
            driver.testField('g23', metric.g23, 3, '2D', stag)
        driver.testField('J', metric.J, 3, '2D', stag)
        driver.testField('G1', metric.G1, 2, '2D', stag)
        driver.testField('G2', metric.G2, 2, '2D', stag)
        driver.testField('G3', metric.G3, 2, '2D', stag)
        exit(0)

    # test operators...
    if args.operator is not None:
        try:
            driver.testOperator(*operator_inputs[args.operator])
        except KeyError:
            try:
                driver.testOperator2(*operator_inputs2[args.operator])
            except KeyError:
                print('Operator '+args.operator+' not found. Available operators for this test are:')
                print(operator_inputs.keys())
                raise
    else:
        for op_input in operator_inputs.values():
            driver.testOperator(*op_input)
        for op_input in operator_inputs2.values():
            driver.testOperator2(*op_input)

    if driver.checkResults():
        print('pass')
        exit(0)
    else:
        print('fail')
        exit(1)
