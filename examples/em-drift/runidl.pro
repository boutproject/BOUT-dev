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

du=file_import("uedge.grd_beta.nc")


data={zeff:fltarr(9), AA:fltarr(9), gam:fltarr(9), omega:fltarr(9), sparn:fltarr(9)}


 i=0
 d=make_d(path="data_1") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar


 i=1
 d=make_d(path="data_2") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=2
 d=make_d(path="data_4") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=3
 d=make_d(path="data_8") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=4
 d=make_d(path="data_16") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=5
 d=make_d(path="data_32") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=6
 d=make_d(path="data_64") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=7
 d=make_d(path="data_128") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 i=8
 d=make_d(path="data_256") & data.zeff[i]=d.zeff & data.AA[i]=d.AA
 res_pproc, d, du, omega=omega, gamma=gamma, manual=manual, spar, wstar, sparsperp=sparsperp, mu=mu
 data.omega[i]=omega & data.gam[i]=gamma & data.sparn[i]=spar/wstar

 getComplexRoots, s, w1 ;-calculate analytic dispersion relation

 d0 = dispersion(mu, sparsperp)
 d1 = dispersion(mu, sparsperp, /e)

 ;STOP

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
   ;image=TVRD(/true)
   ;file='drift_inst_test.jpg'
   ;print, "Saving in file ", file
   ;write_jpeg, file, image,/true


;-compare with original test results (grendel, 5-feb-2007)
   ;RESTORE, 'orig_test2.idl.dat'
   ;error1=max(abs(data.gam-gam_orig)/gam_orig) + max(abs(data.omega-omega_orig)/omega_orig)
   ;print, "Deviation from original test result is" , error1*1e2, " %"


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
exit, status=status
