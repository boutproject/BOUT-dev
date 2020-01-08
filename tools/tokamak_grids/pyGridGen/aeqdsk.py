#!/usr/bin/env python3

from __future__ import print_function
from builtins import str
from builtins import object
import re
import numpy

"""
@brief A-Eqdsk reader class
@version $Id$

Copyright &copy; 2006-2008, Tech-X Corporation, Boulder, CO
See LICENSE file for conditions of use.

The official document describing a-eqdsk files:
http://fusion.gat.com/THEORY/efit/a_eqdsk.html
"""

class Aeqdsk(object):
   
      def __init__(self):
	"""
	Constructor
	"""
     	self.data = {}

      def openFile(self, filename):
	"""
	open aeqdsk file and parse its content
	"""
	
	fmt_1060 = r'^\s*\*\s*([\w\.\-]+)\s+(\d+)\s+(\d+)\s([\w]+)\s+(\d+)\s+(\d+)\s([\w ]+)\s+\d+\s+\d+\s*$'
	fmt_1040 = r'^\s*' + 4*r'([\s\-]\d+\.\d+[Ee][\+\-]\d\d)'
	fmt_1041 = r'^' + 4*r'\s+([ \-]\d+)' 

	lines =  open(filename, 'r').readlines()

	counter = 0
	m = None
	while m == None:
		line = lines[counter]
		m = re.match(fmt_1060, line)
		counter += 1
# read (neqdsk,1060) time(jj),jflag(jj),lflag,limloc(jj), mco2v,mco2r,qmflag
	if m:
		self.data['time'] = float(m.group(1)), 'time ms'
		self.data['jflag'] = int(m.group(2)), '0 if error'
		self.data['lflag'] = int(m.group(3)), '>0 if error'
		self.data['limloc'] = m.group(4), 'IN/OUT/TOP/BOT: limiter inside/outside/top/bot SNT/SNB: single null top/bottom DN: double null'
		self.data['mco2v'] = int(m.group(5)), 'number of vertical CO2 density chords'
		self.data['mco2r'] = int(m.group(6)), 'number of radial CO2 density chords'
		self.data['qmflag'] = m.group(7), 'axial q(0) flag, FIX if constrained and CLC for float'
	else:
	        raise 'Read error at line %d' % (counter-1)

# read (neqdsk,1040) tsaisq(jj),rcencm,bcentr(jj),pasmat(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['tsaisq'] = float(m.group(1)), "total chi2 from magnetic probes, flux loops, Rogowskiand external coils"
	   self.data['rcencm'] = float(m.group(2)), "major radius in cm for vacuum field BCENTR" 
	   self.data['bcentr'] = float(m.group(3)), "vacuum toroidal magnetic field in Tesla at RCENCM"
	   self.data['pasmat'] = float(m.group(4)), "measured plasma toroidal current in Ampere"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) cpasma(jj),rout(jj),zout(jj),aout(jj)  
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['cpasma'] = float(m.group(1)),  "fitted plasma toroidal current in Ampere-turn"
	   self.data['rout'] = float(m.group(2)),  "major radius of geometric center in cm"
	   self.data['zout'] = float(m.group(3)), "Z of geometric center in cm"
	   self.data['aout'] = float(m.group(4)), "plasma minor radius in cm"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)

# read (neqdsk,1040) eout(jj),doutu(jj),doutl(jj),vout(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['eout'] = float(m.group(1)),  "Plasma boundary elongation"
	   self.data['doutu'] = float(m.group(2)),  "upper triangularity"
	   self.data['doutl'] = float(m.group(3)), "lower triangularity"
	   self.data['vout'] = float(m.group(4)), "plasma volume in cm3"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) rcurrt(jj),zcurrt(jj),qsta(jj),betat(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['rcurrt'] = float(m.group(1)),  "major radius in cm of current centroid"
	   self.data['zcurrt'] = float(m.group(2)),  "Z in cm at current centroid"
	   self.data['qsta'] = float(m.group(3)), "equivalent safety factor q*"
	   self.data['betat'] = float(m.group(4)), "toroidal b in %"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) betap(jj),ali(jj),oleft(jj),oright(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['betap'] = float(m.group(1)),  "poloidal b with normalization average poloidal magnetic BPOLAV defined through Ampere's law"
	   self.data['ali'] = float(m.group(2)),  "li with normalization average poloidal magnetic defined through Ampere's law"
	   self.data['oleft'] = float(m.group(3)), "plasma inner gap in cm"
	   self.data['oright'] = float(m.group(4)), "plasma outer gap in cm"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)

	
# read (neqdsk,1040) otop(jj),obott(jj),qpsib(jj),vertn(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['otop'] = float(m.group(1)),  "plasma top gap in cm"
	   self.data['obott'] = float(m.group(2)),  "plasma bottom gap in cm"  
	   self.data['qpsib'] = float(m.group(3)),  "q at 95% of poloidal flux"
	   self.data['vertn'] = float(m.group(4)),  "vacuum field index at current centroid"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)

	mco2v = self.data['mco2v'][0]
	print('mco2v=', mco2v)

# read (neqdsk,1040) (rco2v(k,jj),k=1,mco2v)
	data = []
	while len(data) < mco2v:
	      line = lines[counter]
	      data += eval('[' + re.sub(r'(\d)\s*([\s\-])\s*(\d)', '\\1, \\2\\3', line) + ']')
	      counter += 1
	self.data['rco2v'] = numpy.array(data),   "path length in cm of vertical CO2 density chord"
	
# read (neqdsk,1040) (dco2v(jj,k),k=1,mco2v)
	data = []
	while len(data) < mco2v:
	      line = lines[counter]
	      data += eval('[' + re.sub(r'(\d)\s*([\s\-])\s*(\d)', '\\1, \\2\\3', line) + ']')
	      counter += 1
	self.data['dco2v'] = numpy.array(data),   "line average electron density in cm3 from vertical CO2 chord"
	
	mco2r = self.data['mco2r'][0]
	print('mco2r=', mco2r)

# read (neqdsk,1040) (rco2r(k,jj),k=1,mco2r)
	data = []
	while len(data) < mco2r:
	      line = lines[counter]
	      data += eval('[' + re.sub(r'(\d)\s*([\s\-])\s*(\d)', '\\1, \\2\\3', line) + ']')
	      counter += 1
	self.data['rco2r'] = numpy.array(data),   "path length in cm of radial CO2 density chord"

# read (neqdsk,1040) (dco2r(jj,k),k=1,mco2r)	
	data = []
	while len(data) < mco2r:
	      line = lines[counter]
  	      data += eval('[' + re.sub(r'(\d)\s*([\s\-])\s*(\d)', '\\1, \\2\\3', line) + ']')
              counter += 1
	self.data['dco2r'] = numpy.array(data),   "line average electron density in cm3 from radial CO2 chord"
	
# read (neqdsk,1040) shearb(jj),bpolav(jj),s1(jj),s2(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['shearb'] = float(m.group(1)),   ""
	   self.data['bpolav'] = float(m.group(2)),   "average poloidal magnetic field in Tesla defined through Ampere's law"
	   self.data['s1'] = float(m.group(3)),  "Shafranov boundary line integrals"
	   self.data['s2'] = float(m.group(4)),  "Shafranov boundary line integrals"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) s3(jj),qout(jj),olefs(jj),orighs(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['s3'] = float(m.group(1)),   "Shafranov boundary line integrals"
	   self.data['qout'] = float(m.group(2)),   "q at plasma boundary"
	   self.data['olefs'] = float(m.group(3)),  ""
	   self.data['orighs'] = float(m.group(4)),  "outer gap of external second separatrix in cm"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) otops(jj),sibdry(jj),areao(jj),wplasm(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['otops'] = float(m.group(1)),   "top gap of external second separatrix in cm"
	   self.data['sibdry'] = float(m.group(2)),   ""
	   self.data['areao'] = float(m.group(3)),  "cross sectional area in cm2"
	   self.data['wplasm'] = float(m.group(4)),  ""
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) terror(jj),elongm(jj),qqmagx(jj),cdflux(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['terror'] = float(m.group(1)),   "equilibrium convergence error"
	   self.data['elongm'] = float(m.group(2)),  "elongation at magnetic axis" 
	   self.data['qqmagx'] = float(m.group(3)),  "axial safety factor q(0)"
	   self.data['cdflux'] = float(m.group(4)),  "computed diamagnetic flux in Volt-sec"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) alpha(jj),rttt(jj),psiref(jj),xndnt(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['alpha'] = float(m.group(1)),   "Shafranov boundary line integral parameter"
	   self.data['rttt'] = float(m.group(2)),   "Shafranov boundary line integral parameter"
	   self.data['psiref'] = float(m.group(3)),  "reference poloidal flux in VS/rad"
	   self.data['xndnt'] = float(m.group(4)),  "vertical stability parameter, vacuum field index normalized to critical index value"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) rseps(1,jj),zseps(1,jj),rseps(2,jj),zseps(2,jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['rseps1'] = float(m.group(1)),   "major radius of x point in cm"
	   self.data['zseps1'] = float(m.group(2)),   ""
	   self.data['rseps2'] = float(m.group(3)),  "major radius of x point in cm"
	   self.data['zseps2'] = float(m.group(4)),  ""
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) sepexp(jj),obots(jj),btaxp(jj),btaxv(jj)	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['sepexp'] = float(m.group(1)),   "separatrix radial expansion in cm"
	   self.data['obots'] = float(m.group(2)),   "bottom gap of external second separatrix in cm"
	   self.data['btaxp'] = float(m.group(3)),  "toroidal magnetic field at magnetic axis in Tesla"
	   self.data['btaxv'] = float(m.group(4)),  "vacuum toroidal magnetic field at magnetic axis in Tesla"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) aaq1(jj),aaq2(jj),aaq3(jj),seplim(jj)	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['aaq1'] = float(m.group(1)),   "minor radius of q=1 surface in cm, 100 if not found"
	   self.data['aaq2'] = float(m.group(2)),   "minor radius of q=2 surface in cm, 100 if not found"
	   self.data['aaq3'] = float(m.group(3)),  "minor radius of q=3 surface in cm, 100 if not found"
	   self.data['seplim'] = float(m.group(4)),  "> 0 for minimum gap in cm in divertor configurations, < 0 absolute value for minimum distance to external separatrix in limiter configurations"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040) rmagx(jj),zmagx(jj),simagx(jj),taumhd(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['rmagx'] = float(m.group(1)),   "major radius in cm at magnetic axis"
	   self.data['zmagx'] = float(m.group(2)),   ""
	   self.data['simagx'] = float(m.group(3)),  ""
	   self.data['taumhd'] = float(m.group(4)),  "energy confinement time in ms"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040,err=380) betapd(jj),betatd(jj),wplasmd(jj),diamag(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['betapd'] = float(m.group(1)),   "diamagnetic poloidal b"
	   self.data['betatd'] = float(m.group(2)),   "diamagnetic toroidal b in %"
	   self.data['wplasmd'] = float(m.group(3)),  "diamagnetic plasma stored energy in Joule"
	   self.data['fluxx'] = float(m.group(4)),  "measured diamagnetic flux in Volt-sec"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# read (neqdsk,1040,err=380) vloopt(jj),taudia(jj),qmerci(jj),tavem(jj)
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['vloopt'] = float(m.group(1)),   "measured loop voltage in volt"
	   self.data['taudia'] = float(m.group(2)),   "diamagnetic energy confinement time in ms"
	   self.data['qmerci'] = float(m.group(3)),  "Mercier stability criterion on axial q(0), q(0) > QMERCI for stability"
	   self.data['tavem'] = float(m.group(4)),  "average time in ms for magnetic and MSE data"
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
# ishot > 91000
	line = lines[counter]
	m = re.match(fmt_1041, line)
	if m:
	   self.data['nsilop'] = int(m.group(1)),   ""
	   self.data['magpri'] = int(m.group(2)),   ""
	   self.data['nfcoil'] = int(m.group(3)),  ""
	   self.data['nesum'] = int(m.group(4)),  ""
	   counter += 1
	else:
	   raise 'Read error at line %d:%s' % (counter, line)
	
	nsilop = self.data['nsilop'][0]
	magpri = self.data['magpri'][0]
	print('nsilop=', nsilop, ' magpri=', magpri)
	data = []
	while len(data) < nsilop + magpri:
	      line = lines[counter]
	      data += eval('[' + re.sub(r'(\d)\s*([\s\-])\s*(\d)', '\\1, \\2\\3', line) + ']')
	      counter += 1
	self.data['csilop'] = numpy.array( data[:nsilop] ), "computed flux loop signals in Weber" 
	self.data['cmpr2'] = numpy.array( data[nsilop:] ), ""
	   
#
	data = []
	nfcoil = self.data['nfcoil'][0]
	while len(data) < nfcoil:
	      line = lines[counter]
	      data += eval('[' + re.sub(r'(\d)\s*([\s\-])\s*(\d)', '\\1, \\2\\3', line) + ']')
	      counter += 1
	self.data['ccbrsp'] = numpy.array(data), "computed external coil currents in Ampere"
	
	data = []
	nesum = self.data['nesum'][0]
	while len(data) < nesum:
	      line = lines[counter]
	      data += eval('[' + re.sub(r'(\d)\s*([\s\-])\s*(\d)', '\\1, \\2\\3', line) + ']')
	      counter += 1
	self.data['eccurt'] = numpy.array(data), "measured E-coil current in Ampere"
	
#
	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['pbinj'] = float(m.group(1)),   "neutral beam injection power in Watts"
	   self.data['rvsin'] = float(m.group(2)),   "major radius of vessel inner hit spot in cm"
	   self.data['zvsin'] = float(m.group(3)),  "Z of vessel inner hit spot in cm"
	   self.data['rvsout'] = float(m.group(4)),  "major radius of vessel outer hit spot in cm"
	counter += 1
	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['zvsout'] = float(m.group(1)),   "Z of vessel outer hit spot in cm"
	   self.data['vsurfa'] = float(m.group(2)),   "plasma surface loop voltage in volt, E EQDSK only"
	   self.data['wpdot'] = float(m.group(3)),  "time derivative of plasma stored energy in Watt, E EQDSK only"
	   self.data['wbdot'] = float(m.group(4)),  "time derivative of poloidal magnetic energy in Watt, E EQDSK only"
	counter += 1
	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['slantu'] = float(m.group(1)),   ""
	   self.data['slantl'] = float(m.group(2)),   ""
	   self.data['zuperts'] = float(m.group(3)),  ""
	   self.data['chipre'] = float(m.group(4)),  "total chi2 pressure"
	counter += 1
	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['cjor95'] = float(m.group(1)),   ""
	   self.data['pp95'] = float(m.group(2)),   "normalized P'(y) at 95% normalized poloidal flux"
	   self.data['ssep'] = float(m.group(3)),  ""
	   self.data['yyy2'] = float(m.group(4)),  "Shafranov Y2 current moment"
	   counter += 1
	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['xnnc'] = float(m.group(1)),   ""
	   self.data['cprof'] = float(m.group(2)),   "current profile parametrization parameter"
	   #self.data['oring'] = float(m.group(3)),  ""  (not used)
	   self.data['cjor0'] = float(m.group(4)),  "normalized flux surface average current density at 99% of normalized poloidal flux"
	   counter += 1
	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['fexpan'] = float(m.group(1)),   "flux expansion at x point"
	   self.data['qqmin'] = float(m.group(2)),   "minimum safety factor qmin"
	   self.data['chigamt'] = float(m.group(3)),  "total chi2 MSE"
	   self.data['ssi01'] = float(m.group(4)),  "magnetic shear at 1% of normalized poloidal flux"
	   counter += 1
	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['fexpvs'] = float(m.group(1)),   "flux expansion at outer lower vessel hit spot"
	   self.data['sepnose'] = float(m.group(2)),   "radial distance in cm between x point and external field line at ZNOSE"
	   self.data['ssi95'] = float(m.group(3)),  "magnetic shear at 95% of normalized poloidal flux"
	   self.data['rqqmin'] = float(m.group(4)),  "normalized radius of qmin , square root of normalized volume"
	   counter += 1
	
	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['cjor99'] = float(m.group(1)),   ""
	   self.data['cj1ave'] = float(m.group(2)),   "normalized average current density in plasma outer 5% normalized poloidal flux region"
	   self.data['rmidin'] = float(m.group(3)),  "inner major radius in m at Z=0.0"
	   self.data['rmidout'] = float(m.group(4)),  "outer major radius in m at Z=0.0"
	   counter += 1
	

	line = lines[counter]
	m = re.match(fmt_1040, line)
	if m:
	   self.data['psurfa'] = float(m.group(1)),   "plasma boundary surface area in m2"
	   #self.data[''] = float(m.group(2)),   ""
	   #self.data[''] = float(m.group(3)),  ""
	   #self.data[''] = float(m.group(4)),  ""
	   counter += 1
	

      def getAll(self):
	return self.data

      def getAllVars(self):
        return list(self.data.keys())

      def get(self, varname): 
	return self.data[varname]

################################

def main():
	import sys
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
                  help="g-eqdsk file", default="")
        parser.add_option("-a", "--all", dest="all",
                  help="display all variables", action="store_true",)
        parser.add_option("-v", "--vars", dest="vars", 
                  help="comma separated list of variables (use '-v \"*\"' for all)", default="*")
        #parser.add_option("-p", "--plot", dest="plot",
        #          help="plot all variables", action="store_true",)
        parser.add_option("-i", "--inquire", dest="inquire",
                  help="inquire list of variables", action="store_true",)

        options, args = parser.parse_args()
	if not options.filename:
	   parser.error("MUST provide filename (type -h for list of options)")

	eq = Aeqdsk()
	eq.openFile(options.filename)
	
	if options.inquire:
	   print(eq.getAllVars())

	if options.all:
	   print(eq.getAll())

	vs = eq.getAllVars()
	if options.vars != '*':
	   vs = options.vars.split(',')

	for v in vs:
	    print('%s: %s'% (v, str(eq.get(v))))



	

if __name__ == '__main__': main()

