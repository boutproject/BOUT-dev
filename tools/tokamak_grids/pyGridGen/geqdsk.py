#!/usr/bin/env python3

from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from builtins import object
import re
import numpy

"""
@brief G-Eqdsk reader class
@version $Id$

Copyright &copy; 2006-2008, Tech-X Corporation, Boulder, CO
See LICENSE file for conditions of use.

The official document describing g-eqdsk files:
https://fusion.gat.com/conferences/snowmass/working/mfe/physics/p3/equilibria/g_eqdsk_s.pdf
"""

class Geqdsk(object):

      def __init__(self):
        """
        Constructor
        """
        self.data = {}

      def openFile(self, filename):
        """
        open geqdsk file and parse its content
        """

        lines =  open(filename, 'r').readlines()

        # first line
        m = re.search(r'^\s*(.*)\s+\d+\s+(\d+)\s+(\d+)\s*$', lines[0])
        self.data['case'] = m.group(1), "Identification character string"
        self.data['nw'] = int(m.group(2)), "Number of horizontal R grid points"
        self.data['nh'] = int(m.group(3)), "Number of vertical Z grid points"

        fltsPat = r'^\s*([ \-]\d\.\d+[Ee][\+\-]\d\d)([ \-]\d\.\d+[Ee][\+\-]\d\d)([ \-]\d\.\d+[Ee][\+\-]\d\d)([ \-]\d\.\d+[Ee][\+\-]\d\d)([ \-]\d\.\d+[Ee][\+\-]\d\d)\s*$'

        # 2nd line
        m = re.search(fltsPat, lines[1])
        self.data['rdim'] = float(m.group(1)), "Horizontal dimension in meter of computational box"
        self.data['zdim'] = float(m.group(2)), "Vertical dimension in meter of computational box"
        self.data['rcentr'] = float(m.group(3)), "R in meter of vacuum toroidal magnetic field BCENTR"
        self.data['rleft'] = float(m.group(4)), "Minimum R in meter of rectangular computational box"
        self.data['zmid'] = float(m.group(5)), "Z of center of computational box in meter"

        # 3rd line
        m = re.search(fltsPat, lines[2])
        self.data['rmaxis'] = float(m.group(1)), "R of magnetic axis in meter"
        self.data['zmaxis'] = float(m.group(2)), "Z of magnetic axis in meter"
        self.data['simag'] = float(m.group(3)), "poloidal flux at magnetic axis in Weber /rad"
        self.data['sibry'] = float(m.group(4)), "poloidal flux at the plasma boundary in Weber /rad"
        self.data['bcentr'] = float(m.group(5)), "Vacuum toroidal magnetic field in Tesla at RCENTR"

        # 4th line
        m = re.search(fltsPat, lines[3])
        self.data['current'] = float(m.group(1)), "Plasma current in Ampere"
        #self.data['simag'] = float(m.group(2)), ""

        #self.data['rmaxis'] = float(m.group(4)), ""

        # 5th line
        m = re.search(fltsPat, lines[4])
        #self.data['zmaxis'] = float(m.group(1)), ""

        #self.data['sibry'] = float(m.group(3)), ""

        # read remaining data
        data = []
        counter = 5
        while 1:
                line = lines[counter]
                m = re.match(r'^\s*[ \-]\d\.\d+[Ee][\+\-]\d\d', line)
                if not m: break
                data += eval('[' + re.sub(r'(\d)([ \-]\d\.)', '\\1,\\2', line) + ']')
                counter += 1

        nw = self.data['nw'][0]
        nh = self.data['nh'][0]

        self.data['fpol'] = numpy.array(data[0:nw]), "Poloidal current function in m-T, F = RBT on flux grid"
        self.data['pres'] = numpy.array(data[nw:2*nw]), "Plasma pressure in nt / m 2 on uniform flux grid"

        self.data['ffprime'] = numpy.array(data[2*nw:3*nw]), "FF'(psi) in (mT)^2/(Weber/rad) on uniform flux grid"

        self.data['pprime'] = numpy.array(data[3*nw:4*nw]), "P'(psi) in (nt/m2)/(Weber/rad) on uniform flux grid"

        self.data['psirz'] = numpy.reshape( data[4*nw:4*nw+nw*nh], (nh, nw) ), "Poloidal flux in Weber / rad on the rectangular grid points"
        self.data['qpsi']  = numpy.array(data[4*nw+nw*nh:5*nw+nw*nh]), "q values on uniform flux grid from axis to boundary"

        line = lines[counter]
        m = re.search(r'^\s*(\d+)\s+(\d+)', line)
        print(line)
        nbbbs = int(m.group(1))
        limitr = int(m.group(2))
        self.data['nbbbs'] = nbbbs, "Number of boundary points"
        self.data['limitr'] = limitr, "Number of limiter points"
        counter += 1

        data = []
        while 1:
                line = lines[counter]
                m = re.search(r'^\s*[ \-]\d\.\d+[Ee][\+\-]\d\d', line)
                counter += 1
                if not m: break
                data += eval('[' + re.sub(r'(\d)([ \-]\d\.)', '\\1,\\2', line) + ']')
        self.data['rbbbs'] = numpy.zeros( (nbbbs,), numpy.float64 ), "R of boundary points in meter"
        self.data['zbbbs'] = numpy.zeros( (nbbbs,), numpy.float64 ), "Z of boundary points in meter"
        for i in range(nbbbs):
                self.data['rbbbs'][0][i] = data[2*i]
                self.data['zbbbs'][0][i] = data[2*i + 1]

        self.data['rlim'] = numpy.zeros( (limitr,), numpy.float64 ), "R of surrounding limiter contour in meter"
        self.data['zlim'] = numpy.zeros( (limitr,), numpy.float64 ), "Z of surrounding limiter contour in meter"
        for i in range(limitr):
                self.data['rlim'][0][i] = data[2*nbbbs + 2*i]
                self.data['zlim'][0][i] = data[2*nbbbs + 2*i + 1]


      def getAll(self):
        return self.data

      def getAllVars(self):
        return list(self.data.keys())

      def get(self, varname):
        return self.data[varname.lower()][0]

      def getDescriptor(self, varname):
        return self.data[varname.lower()][1]

################################

def main():
        from optparse import OptionParser
        parser = OptionParser()
        parser.add_option("-f", "--file", dest="filename",
                  help="g-eqdsk file", default="")
        parser.add_option("-a", "--all", dest="all",
                  help="display all variables", action="store_true",)
        parser.add_option("-v", "--vars", dest="vars",
                  help="comma separated list of variables (use '-v \"*\"' for all)", default="*")
        parser.add_option("-p", "--plot", dest="plot",
                  help="plot all variables", action="store_true",)
        parser.add_option("-i", "--inquire", dest="inquire",
                  help="inquire list of variables", action="store_true",)


        options, args = parser.parse_args()
        if not options.filename:
           parser.error("MUST provide filename (type -h for list of options)")


        geq = Geqdsk()
        geq.openFile(options.filename)

        if options.inquire:
           print(geq.getAllVars())

        if options.all:
           print(geq.getAll())

        vs = geq.getAllVars()
        if options.vars != '*':
           vs = options.vars.split(',')

        for v in vs:
            print('%s: %s'% (v, str(geq.get(v))))

        if options.plot:
           from matplotlib import pylab

           if options.vars == '*':
              options.vars = geq.getAllVars()
              print(options.vars)
           else:
              vs = options.vars.split(',')
              options.vars = vs

           xmin = geq.get('simag')
           xmax = geq.get('sibry')
           nx   = geq.get('nw')
           dx = float(xmax - xmin) / float(nx - 1)
           x = numpy.arange(xmin, xmin + (xmax-xmin)*(1.+1.e-6), dx)
           for v in options.vars:
               if v[0] != 'r' and v[0] != 'z':
                  data = geq.get(v)
                  if len(numpy.shape(data)) == 1:
                     pylab.figure()
                     pylab.plot(x, data)
                     pylab.xlabel('psi poloidal')
                     pylab.ylabel(v)
                     pylab.title(geq.getDescriptor(v))
           # 2d plasma plot
           nw = geq.get('nw')
           nh = geq.get('nh')
           rmin = geq.get('rleft')
           rmax = rmin + geq.get('rdim')
           dr = float(rmax - rmin)/float(nw - 1)
           zmin = geq.get('zmid') - geq.get('zdim')/2.0
           zmax = geq.get('zmid') + geq.get('zdim')/2.0
           dz = (zmax - zmin)/float(nh - 1)
           rs = numpy.arange(rmin, rmin + (rmax-rmin)*(1.+1.e-10), dr)
           zs = numpy.arange(zmin, zmin + (zmax-zmin)*(1.+1.e-10), dz)
           pylab.figure()
           pylab.pcolor(rs, zs, geq.get('psirz'), shading='interp')
           pylab.plot(geq.get('rbbbs'), geq.get('zbbbs'), 'w-')
           #pylab.plot(geq.get('rlim'), geq.get('zlim'), 'k--')
           pylab.axis('image')
           pylab.title('poloidal flux')
           pylab.xlabel('R')
           pylab.ylabel('Z')

           pylab.show()


if __name__ == '__main__': main()
