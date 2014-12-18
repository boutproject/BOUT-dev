import numpy as np
from boutdata import collect
from boututils import moment_xyzt, file_import
from pylab import save, figure, plot, title, xlabel, ylabel
from elm_size import elm_size

path='/Users/brey/BOUT/Nonlinear_n15_lund1e5.George/data'

t_array=collect('t_array', path=path)
save('t_array.dat', t_array)
p0=collect('P0', path=path)
save('p0.dat', p0)


# n0=collect('n0', path=path)
# save('n0.dat', n0
# ti0=collect('ti0', path=path)
# save('ti0.dat', ti0)
# te0=collect('te0', path=path)
# save('te0.idl.dat', te0)

gfile=file_import('/Users/brey/BOUT/cbm18_6_y064_x516_102709.grd.cdl')

p=collect('P', path=path,tind=[0,1240])
save('p.dat', p)
res=moment_xyzt(p,'RMS','DC')
rmsp=res.rms
dcp=res.dc
save('rmsp.dat', rmsp)
save('dcp.dat',  dcp)
# gfile=file_import('data/cbm18_dens6.x516_y64.nc')
elmsp=elm_size(dcp,p0,gfile,yind=32,Bbar=gfile['bmag'])
save('elmsp.dat',  elmsp)

figure(0)
plot(t_array,elmsp.s2, 'k-')
xlabel('t/Ta')
ylabel('Elm size')
title('Elm size, P')


phi=collect('phi', path=path ,tind=[0,1240])
save('phi.dat', phi)
res=moment_xyzt( phi, 'DC', 'RMS')
save('dcphi.dat',res.dc)
save('rmsphi.dat', res.rms)

