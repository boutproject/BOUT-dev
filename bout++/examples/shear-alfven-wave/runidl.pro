;pathnow=getenv('PWD')
;cd,getenv('BOUT_DA')
;.run PDB2IDL/pdb2idl.pro
;.run Plots2D/moment_xyzt.pro
;.run Plots2D/allrms.pro
;cd,pathnow

device, decomposed=0
tek_color

.r pdb2idl
.r moment_xyzt
.r allrms

;safe_colors, /first

 fu='uedge.grd_Te_10.nc'
 
 path = "data_0.15/"
 @alf_pproc
 k015=kPerp & v015=vPulse

 path = "data_0.25/"
 @alf_pproc
 k025=kPerp & v025=vPulse

 path = "data_0.50/"
 @alf_pproc
 k050=kPerp & v050=vPulse

 path = "data_1.00/"
 @alf_pproc
 k100=kPerp & v100=vPulse

 vte=0.0
 va=2.2e8  ;cm/s @1 T, 1e14 cm-3
 wpe=5.6e11  ;rad/s @1 T
 c=3e10    ;cm/s

 safe_colors, /first

 plot, [k015,k025,k050,k100], [v015,v025,v050,v100], psym=4, yr=[0,3e8], syms=3,$
   xtit='kperp, cm-1', ytit='V, cm/s', tit='Dispersion relation for Alfven wave', col=1, chars=2
 kk=35*findgen(101)/100.
 oplot, kk, va*sqrt((1+((vte/va)*(kk*c/wpe))^2)/(1+(kk*c/wpe)^2)), col=1, lin=2

 ;-save the picture
 ;image=TVRD(/true)
 ;file='alftest_noPe.jpg'
 ;print, "Saving in file ", file
 ;write_jpeg, file, image,/true

 ;-compare results with original test
 restore, 'orig_test.idl.dat'
 result_31jan2007=[v015orig, v025orig, v050orig, v100orig]
 error1=max(abs(result_31jan2007-[v015,v025,v050,v100])/result_31jan2007) 
 print, "Deviation from original test result is" , error1*1e2, " %"

 oplot, [k015,k025,k050,k100], [v015orig, v025orig, v050orig, v100orig], col=1, psym=7


 ;-compare results with analytic answer
 result_analytic=[1.08129e+08,  1.50732e+08,  1.94260e+08,  2.12618e+08]
 error2=max(abs(result_analytic-[v015,v025,v050,v100])/result_analytic) 
 print, "Deviation from analytic result is" , error2*1e2, " %"


 ;- Output to a Postscript file
 SET_PLOT, 'PS'
 DEVICE, file="shearwave.ps"
 plot, [k015,k025,k050,k100], [v015,v025,v050,v100], psym=4, yr=[0,3e8], syms=3,$
   xtit='kperp, cm-1', ytit='V, cm/s', tit='Dispersion relation for Alfven wave', col=1
 oplot, kk, va*sqrt((1+((vte/va)*(kk*c/wpe))^2)/(1+(kk*c/wpe)^2)), col=1, lin=2
 oplot, [k015,k025,k050,k100], [v015orig, v025orig, v050orig, v100orig], col=1, psym=7
 DEVICE, /close
 SET_PLOT, 'X'

 for i=10,0,-1 do begin print, "exiting in ", i, " seconds" & wait,1
 
 if (error1 gt 1e-4) then status=1 else status=0
 print, 'status=', status
 exit, status=status
