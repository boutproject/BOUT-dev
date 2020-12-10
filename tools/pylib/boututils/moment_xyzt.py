from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np
from boututils.bunch import Bunch

def RMSvalue( vec1d):
#;
#; -get rms of a 1D signal
#;------------------------

    nel=np.size(vec1d)
    valav=old_div(np.sum(vec1d),nel) 
    valrms=np.sqrt(old_div(np.sum((vec1d-valav)**2),nel))
    acvec=vec1d-valav

    return Bunch(valrms=valrms,
                 valav=valav,
                 acvec=acvec)



def moment_xyzt( sig_xyzt, *args):#rms=None, dc=None, ac=None):
#;
#; Calculate moments of a 4d signal of (x,y,z,t), i.e,
#; -RMS, i.e., a function of (x,y,t)
#; -DC (average in z), i.e., a function of (x,y,t)
#; -AC (DC subtracted out), i.e., a function of (x,y,z,t)
#;-------------------------------------------------------------------

        d = np.shape(sig_xyzt)
        if np.size(d) != 4 :
            print("Error: Variable must be 4D (x,y,z,t)")
            return
            
    
        siz=np.shape(sig_xyzt)
        rms=np.zeros((siz[0],siz[1],siz[2]))
        dc=np.zeros((siz[0],siz[1],siz[2]))
        if 'AC' in args : ac=np.zeros((siz[0],siz[1],siz[2],siz[3]))
  
        
        data = sig_xyzt
        if np.modf(np.log2(siz[3]))[0] != 0.0 :
            print("WARNING: Expecting a power of 2 in Z direction")
      
            if np.modf(np.log2(siz[3]-1))[0] and (siz[3] > 1) :
                print(" -> Truncating last point to get power of 2")
                data = data[:,:,0:(siz[3]-2),:]
                siz[3] = siz[3] - 1
              
  
        for ix in range (siz[1]):
            for iy in range (siz[2]):
                for it in range (siz[0]):
                    val=RMSvalue(sig_xyzt[it,ix,iy,:])
              
                    rms[it,ix,iy]=val.valrms
                    dc[it,ix,iy]=val.valav
                    if 'AC' in args : ac[it,ix,iy,:]=[val.acvec,val.acvec[0]]
         
        res=Bunch()
         
        if 'RMS' in args:
            res.rms = rms
        if 'DC' in args:
            res.dc = dc
        if 'AC' in args:
            res.ac = ac

        if 'RMS' not in args and 'DC' not in args and 'AC' not in args :
            raise RuntimeError('Wrong argument')
        return res
