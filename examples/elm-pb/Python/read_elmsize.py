from boutdata.collect import collect
from boututils.moment_xyzt import moment_xyzt
from boututils.datafile import DataFile
from pylab import save, figure, plot, title, xlabel, ylabel, show, tight_layout
from elm_size import elm_size

path='./data'

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

with DataFile("./cbm18_dens8.grid_nx68ny64.nc") as f:
    gfile = {v: f.read(v) for v in f.keys()}

p=collect('P', path=path)
save('p.dat', p)
res=moment_xyzt(p,'RMS','DC')
rmsp=res.rms
dcp=res.dc
save('rmsp.dat', rmsp)
save('dcp.dat',  dcp)
elmsp=elm_size(dcp,p0,gfile,yind=32,Bbar=gfile['bmag'])
save('elmsp.dat',  elmsp)

figure(0)
plot(t_array,elmsp.s2, 'k-')
xlabel('t/Ta')
ylabel('Elm size')
title('Elm size, P')
tight_layout()
show()


phi=collect('phi', path=path )
save('phi.dat', phi)
res=moment_xyzt( phi, 'DC', 'RMS')
save('dcphi.dat',res.dc)
save('rmsphi.dat', res.rms)
