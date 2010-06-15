pro res_pproc, d, du, omega=omega, gamma=gamma, spar, wstar, MANUAL=MANUAL, DEBUG=DEBUG, skip=skip, $
               sparsperp=sparsperp, mu=mu
;
; Postprocessing of resisitive instability test data
;
;-------------------------------------------------------------;

IF NOT KEYWORD_SET(skip) THEN skip = 15 ; skip first part of the curve

;tek_color
safe_colors, /first

 nstep=1
 nt=n_elements(d.ni_xyzt[0,0,0,*])
 nx=n_elements(d.ni_xyzt[*,0,0,0])
 ny=n_elements(d.ni_xyzt[0,*,0,0])
 nz=n_elements(d.ni_xyzt[0,0,*,0])


 ;;-first calculate geometric and physical quantities

      lZeta=d.zmax*(d.rho_s*1e2)*2*!PI*(du.R0/du.hthe0) ;-toroidal range [cm]
      lbNorm=lZeta*(du.BPXY[0,ny/2]/du.BXY[0,ny/2])    ;-binormal coord range [cm]
      zPerp=lbNorm*findgen(nz)/(nz-1) ;-binormal coordinate [cm]

      cLight=3e10 ;-speed of light [cm/s]
      vTe=4.2e7*sqrt(du.Te_x)   ;-electron thermal speed [cm/s]
      kperp=2*!PI/max(zPerp)    ;-binormal wavenumber, [cm-1]
      wce=1.76e7*1e4*du.Bmag    ;-electron cyclotron frequency, [rad/s]

      Ln=MEAN(ABS(du.ni0[*,ny/2]/DERIV(du.Rxy[*,ny/2]*1e2,du.ni0[*,ny/2]))) ;-Ni scale length [cm]
      vPe=(vTe)^2/(wce*Ln) ;-electron diamagnetic drift speed [cm/s]
      wstar=vPe*kPerp  
      print, "wstar=", wstar, " [rad/s]"

      logLam=24. - alog(sqrt(du.ni_x*1e14)/du.te_x)
      nuei=d.zeff*2.91e-6*(du.ni_x*1e14)*logLam/(du.te_x)^1.5
      wci=9.58e3*(1./d.AA)*1e4*du.Bmag ;-ion cyclotron frequency for Mi/mp=d.AA, [rad/s]

      lpar=total((du.bxy/du.bpxy)*du.dlthe)/du.nx ;-[m], average over flux surfaces
      kpar=2*!pi/(1e2*lpar) ;cm-1
      spar=(kpar/kperp)^2 * wci * wce / (0.51 * nuei) ;[1/s]


      wpe=5.64e4*sqrt(1e14*du.ni_x) ;-electron plasma frequency, [rad/s]
      mu=(cLight*kperp/wpe)^2 
      sperp=(0.51 * nuei) * mu ;[1/s]
      print, 'sPerp/wstar=', sperp/wstar
      PRINT, "sPar /wstar=", spar / wstar

      sparsperp = spar*sperp / (wstar*wstar)
      PRINT, "spar * sperp / wstar^2= ", sparsperp

      print, 'mu=', mu
 ;;------------------------------------------------------------------;;


!p.title='Resistive drift instability in BOUT'
!p.title=!p.title+', Zeff=' + strtrim(string(d.zeff,f='(f5.1)'),2)

yscale=2*max([max(d.phi_xyzt[*,20,*,*]), max(d.ni_xyzt[*,20,*,*])])

!p.multi=[0,1,2,0,0]
  xtit='binormal coord'
  ytit='phi'
   plot, d.phi_xyzt[0,20,*,nt-1],/xst, xtit=xtit, ytit=ytit, yr=[-1,1]*yscale,/yst,/nod, chars=1.5
   for i=0,nt-1 do oplot, d.phi_xyzt[0,20,*,i], col=1+(i mod 4)
   oplot, d.phi_xyzt[0,20,*,nt-1]*0

  xtit='binormal coord'
  ytit='ni'
   plot, d.ni_xyzt[0,20,*,nt-1],/xst, xtit=xtit, ytit=ytit, yr=[-1,1]*yscale,/yst,/nod, chars=1.5
   for i=0,nt-1,nStep do oplot, d.ni_xyzt[0,20,*,i], col=1+(i mod 4)
   oplot, d.ni_xyzt[0,20,*,nt-1]*0
!p.multi=0



WAIT,3
 nt0=skip ;-skip initial part of the curve
 maxVal=fltarr(nt-nt0)
 for i=nt0,nt-1 do maxVal[i-nt0]=max(d.ni_xyzt[0,20,*,i])
 plot, d.t_array[nt0:*]/d.wci, alog(maxVal), psym=4, syms=3, xtit='time, s', ytit='ln<Ni>',/yst, chars=1.5

 if keyword_set(MANUAL) then begin
   print, "Mark 2 points on straight line to measure the exponential growth rate"
   print, "Click point 1" & mark, x=x1, y=y1 
   print, "Click point 2" & mark,x=x2,y=y2 
   oplot, [x1,x2], [y1,y2], col=2
   gamma=(y2-y1)/(x2-x1)
 endif else begin
   xx=d.t_array[nt0:*]/d.wci & yy=ALOG(maxVal)
   res=Poly_Fit(xx,yy,1)
   oplot, xx, res[0] + res[1]*xx, col=2
   gamma=res[1]
 endelse

 print, "gamma=", gamma, " s-1"
 ;STOP

WAIT,3 
 xtit="Binormal coordinate, cm"
 ytit="Normalized perturbation"

 plot, zPerp, d.ni_xyzt[0,20,*,nt-1]/maxVal[nt-1-nt0],/xst, yr=[-2,2], xtit=xtit, ytit=ytit, chars=1.5
 for i=nt0,nt-1,nStep do begin 
  oplot, zPerp, d.ni_xyzt[0,20,*,i]/maxVal[i-nt0], col=1+(i mod 4) 
  wait,0.1
 endfor


WAIT,3
;-track the motion of the peak to infer phase velocity

 xPeak=fltarr(nt-nt0)

 for i=nt0,nt-1 do begin 
  sig=d.ni_xyzt[0,20,*,i]/maxVal[i-nt0]  
  zPerp2=min(zPerp) + (max(zPerp)-min(zPerp))*findgen(1001)/1000.
  sig2=spline(zPerp,sig,zPerp2)

  iMax=where(sig2 eq max(sig2)) 
  xPeak[i-nt0]=zPerp2[iMax[0]]
 endfor


 ;-truncate the set to avoid periodic recurrence
 if (xPeak[1] lt xPEak[0]) then begin
   iValid=where(xPeak ge min(zPerp) and xPeak le xPeak[0])
 endif else begin
   iValid=where(xPeak le max(zPerp) and xPeak ge xPeak[0])
 endelse

 xPeak=xPeak(iValid)
 tt=d.t_array(iValid)/d.wci

 plot, tt, xPeak, psym=4, syms=3, chars=1.5, xtit='time, s', ytit='Binormal coordinate of wave peak'


 if keyword_set(MANUAL) then begin
   print, "Marking 2 points on the line to measure the phase velocity"
   print, "Click point 1" & mark, x=x1, y=y1 
   print, "Click point 2" & mark,x=x2,y=y2 
    Vphase=abs((y2-y1)/(x2-x1))
    oplot, [x1,x2], [y1,y2], col=3
 endif else begin
   xx=tt & yy=xPeak
   res=Poly_Fit(xx,yy,1)
   oplot, xx, res[0] + res[1]*xx, col=2
   Vphase=abs(res[1])
 endelse

 print, "Vphase=", Vphase, " cm/s"


 ;-calculate normalized quantities-
  omega=Vphase*kperp/wstar ;-normalized real frequency
  gamma=gamma/wstar ;-normalized growth rate
  print, "Normalized omega=", omega
  print, "Normalized gamma=", gamma
  wait,3

!p.title=''
if keyword_set(DEBUG) then STOP
end
