from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np

def elm_size(dcp,p0,uedge,xmin=None,xmax=None,yind=None,Bbar=None):
  
    lis=[dcp,p0,uedge]
    if np.size(lis) != 3 :
        print("lack of parameters")
        return 0
    

    if xmin == None : xmin=0
    if xmax == None : xmax=327 
    if yind == None : yind=63 # choose the poloidal location for 1D size
    if Bbar == None : Bbar=1.992782  # the normalized magnetic field 

    mydcp=dcp
    myp0=p0
    g=uedge

    PI = 3.1415926
    MU0 = 4.0e-7*PI

    s=np.shape(mydcp)

    if np.ndim(mydcp) != 3 :
        print("dcp should be 3D(t,x,y)")
 

    nt=s[0]
    nx=s[1]
    ny=s[2]

    Dtheta=g['dy']     #using correct poloidal angle
    psixy=g['psixy']
    R=g['Rxy']
    Bp=g['Bpxy']
    hthe=g['hthe']

    Dpsi=np.zeros((nx,ny))
    Dpsi[0,:]=psixy[1,:]-psixy[0,:]
    Dpsi[nx-1,:]=psixy[nx-1,:]-psixy[nx-2,:]
    for i in range(1,nx-2):
        Dpsi[i,:]=old_div((psixy[i+1,:]-psixy[i-1,:]),2)
     

    Ddcp1=np.zeros(nt)
    Ddcp2=np.zeros(nt)
    Ddcp3=np.zeros(nt)
    Tp01=0.
    Tp02=0.
    Tp03=0.

    for t in range(nt) :
        Ddcp3[t]=2.0*PI*np.sum(mydcp[t,xmin:xmax,:]*hthe[xmin:xmax,:]*Dtheta[xmin:xmax,:]*Dpsi[xmin:xmax,:]/Bp[xmin:xmax,:])
        Ddcp2[t]=np.sum(mydcp[t,xmin:xmax,:]*hthe[xmin:xmax,:]*Dtheta[xmin:xmax,:]*Dpsi[xmin:xmax,:]/(R[xmin:xmax,:]*Bp[xmin:xmax,:]))
        Ddcp1[t]=np.sum(mydcp[t,xmin:xmax,yind]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind])) 
    

    Tp03=2.0*PI*np.sum(myp0[xmin:xmax,:]*hthe[xmin:xmax,:]*Dtheta[xmin:xmax,:]*Dpsi[xmin:xmax,:]/Bp[xmin:xmax,:])
    Tp02=np.sum(myp0[xmin:xmax,:]*hthe[xmin:xmax,:]*Dtheta[xmin:xmax,:]*Dpsi[xmin:xmax,:]/(R[xmin:xmax,:]*Bp[xmin:xmax,:]))
    Tp01=np.sum(myp0[xmin:xmax,yind]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind]))

    s1=np.zeros(nt)
    s2=np.zeros(nt)
    s3=np.zeros(nt)
    E_loss=np.zeros(nt)

    s1=old_div(-Ddcp1,Tp01)   #1D elm size
    s2=old_div(-Ddcp2,Tp02)   #2D elm size
    s3=old_div(-Ddcp3,Tp03)   #3D elm size

    E_loss=-Ddcp3*(0.5*Bbar*Bbar/MU0)    #energy loss, unit J
    E_total=Tp03*(0.5*Bbar*Bbar/MU0)     #total energy, unit J

    elmsize=Bunch()
    elmsize.s1=s1
    elmsize.s2=s2
    elmsize.s3=s3
    elmsize.E_loss=E_loss
    elmsize.E_total=E_total

    return elmsize
  
    
