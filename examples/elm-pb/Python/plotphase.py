from __future__ import division
from past.utils import old_div

from numpy import save, load, angle
import matplotlib.pyplot as plt

fphi = load('fphi.npy')

fte = load('fte.npy')
phase_te = angle(old_div(fphi,fte))
save('phase_te', phase_te)
plt.figure()
plt.plot(mean(mean(phase_te[:,:,3,:],axis=1),axis=1))
plt.title('Te')
plt.savefig('image/phase_te.png')
plt.savefig('image/phase_te.eps')

fti = load('fti.npy')
phase_ti = angle(old_div(fphi,fti))
save('phase_ti', phase_ti)
plt.figure()
plt.plot(mean(mean(phase_ti[:,:,3,:],axis=1),axis=1))
plt.title('ti')
plt.savefig('image/phase_ti.png')
plt.savefig('image/phase_ti.eps')

fni = load('fni.npy')
phase_ni = angle(old_div(fphi,fni))
save('phase_ni', phase_ni)
plt.figure()
plt.plot(mean(mean(phase_ni[:,:,3,:],axis=1),axis=1))
plt.title('ni')
plt.savefig('image/phase_ni.png')
plt.savefig('image/phase_ni.eps')

plt.show()

