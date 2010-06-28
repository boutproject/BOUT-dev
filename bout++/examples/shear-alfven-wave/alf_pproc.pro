
du=file_import(fu)

wci = collect(var="wci", path=path)
rho_s = collect(var="rho_s", path=path)
zmax = collect(var="zmax", path=path)
tt = collect(var="t_array", path=path) / wci
jpar = collect(var="jpar", path=path, x=2, z=10)


tt=REFORM(tt)
sig=REFORM(jpar[0,32,0,*])
xpar=10.*(0.5-findgen(64)/63.)

tt2=max(tt)*findgen(1000)/999.
sig2=spline(tt,sig,tt2)

!p.multi=[0,1,2,0,0]
  ytit='poloidal coordinate, m'
  tit='Pulse bouncing'
  contour, TRANSPOSE(reform(jpar[0,*,0,*])),tt,xpar,nlev=11,$
  /xst,/yst,/fil, xtit=xtit, ytit=ytit, tit=tit
  ;
  xtit='time, s'
  tit='Jpar @ x=x0'
  plot, tt, sig,/xst,psym=3, xtit=xtit, tit=tit
  oplot, tt2, sig2,col=3
  oplot, tt, sig,psym=3
  ;
  ;print, 'Click on 1st peak:'
  ;mark, x=t0, y=y0
  ;makesym,s=1
  ;for i=1,10 do plots, i*t0, y0, psym=8, syms=2, col=2, thick=3
  Find_Peaks, sig2, tt2, tMax=t0
!p.multi=0


Bpol=0.5*(min(du.bpxy)+max(du.bpxy))
Btor=0.5*(min(du.btxy)+max(du.btxy))
height=max(du.zxy)-min(du.zxy)
Lpar=height*(Btor/Bpol)*1e2 ;[cm]
vPulse=Lpar/t0
print, 'vPulse = ', vPulse, " cm/s", format='(A10,e12.6,A5)'

R0=du.R0
hthe0=du.hthe0
rho_s=rho_s*1e2 ;[cm]
Ltor=Zmax*2*!pi*R0*100. ;*rho_s/hthe0
Lperp=Ltor*Bpol/Btor
kperp=2*!Pi/Lperp
print, 'kperp=', kperp, ' cm-1'

WAIT, 3
