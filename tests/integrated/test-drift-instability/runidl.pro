;pathnow=getenv('PWD')
;cd,getenv('BOUT_DA')
;.run PDB2IDL/pdb2idl.pro
;.run Plots2D/moment_xyzt.pro
;.run Plots2D/allrms.pro
;cd,pathnow

.run res_pproc
.run complexRoots
.run make_d


WINDOW, 0, XSIZE=500, YSIZE=800

du=file_import("uedge.grd_std.cdl")


data={zeff:fltarr(9), AA:fltarr(9), gam:fltarr(9), omega:fltarr(9), sparn:fltarr(9)}


 i=0 & data.zeff[i]=1. & data.AA[i]=2e0
 d=make_d(path="data_1") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=1 & data.zeff[i]=2. & data.AA[i]=2e0
 d=make_d(path="data_2") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=2 & data.zeff[i]=4. & data.AA[i]=2e0
 d=make_d(path="data_4") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=3 & data.zeff[i]=8. & data.AA[i]=2e0
 d=make_d(path="data_8") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=4 & data.zeff[i]=16. & data.AA[i]=2e0
 d=make_d(path="data_16") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=5 & data.zeff[i]=32. & data.AA[i]=2e0
 d=make_d(path="data_32") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=6 & data.zeff[i]=64. & data.AA[i]=2e0
 d=make_d(path="data_64") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=7 & data.zeff[i]=128. & data.AA[i]=2e0
 d=make_d(path="data_128") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=8 & data.zeff[i]=256.  & data.AA[i]=2e0
 d=make_d(path="data_256") & d.zeff=data.zeff[i] & d.AA=data.AA[i]
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar


 getComplexRoots, s, w1 ;-calculate analytic dispersion relation

 !p.multi=[0,1,2,0,0]
    tit="Resistive drift instability in BOUT"
    xtit="!7r!3!I||!N/!7x!3!I*!N"	 
    ytit="Im(!7x/x!I*!N!3)"
    plot, data.sparn, data.gam,/xl, psym=4, syms=3, yr=[0,0.5], $
     xtit=xtit, ytit=ytit, tit=tit, chars=1.5
    oplot, s, imaginary(w1), lin=2, col=3

    tit="Resistive drift instability in BOUT"
    xtit="!7r!3!I||!N/!7x!3!I*!N"	 
    ytit="Re(!7x/x!I*!N!3)"
    plot, data.sparn, data.omega,/xl, psym=4, syms=3, yr=[0,1.5], $
     xtit=xtit, ytit=ytit, tit=tit, chars=1.5
    oplot, s, float(w1), lin=2, col=3
 !p.multi=0


;-save the picture
   image=TVRD(/true)
   file='drift_inst_test.jpg'
   print, "Saving in file ", file
   write_jpeg, file, image,/true


;-compare with original test results (grendel, 5-feb-2007)
   RESTORE, 'orig_test2.idl.dat'
   error1=max(abs(data.gam-gam_orig)/gam_orig) + max(abs(data.omega-omega_orig)/omega_orig)
   print, "Deviation from original test result is" , error1*1e2, " %"


;-compare results with analytic answer
   arg=data.sparn
   gval=INTERPOL(imaginary(w1),s,arg)
   oval=INTERPOL(float(w1),s,arg)
   error2=MAX(abs(gval-data.gam)/gval) + MAX(abs(oval-data.omega)/oval)
   print, "Deviation from analytic result is" , error2*1e2, " %"



 for i=10,0,-1 do begin print, "exiting in ", i, " seconds" & wait,1
 
 if (error1 gt 1e-4) then status=1 else status=0
 ;status=0
 print, 'status=', status

STOP

;exit, status=status
