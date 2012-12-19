FUNCTION elm_size,dcp,p0,uedge,xmin=xmin,xmax=xmax,yind=yind,Bbar=Bbar
  
Nparameter=N_Params()
IF Nparameter NE 3 THEN BEGIN
  PRINT,"lack of parameters"
  RETURN, 0
ENDIF

IF NOT KEYWORD_SET(xmin) THEN xmin=0
IF NOT KEYWORD_SET(xmax) THEN xmax=327 
IF NOT KEYWORD_SET(yind) THEN yind=64 ; choose the poloidal location for 1D size
IF NOT KEYWORD_SET(Bbar) THEN Bbar=1.992782  ; the normalized magnetic field 

mydcp=dcp
myp0=p0
g=uedge

PI = 3.1415926
MU0 = 4.0e-7*PI

s=size(mydcp)

IF s[0] NE 3 THEN BEGIN
  PRINT,"dcp should be 3D(x,y,t)"
ENDIF

nx=s[1]
ny=s[2]
nt=s[3]

Dtheta=g.dy     ;using correct poloidal angle
psixy=g.psixy
R=g.rxy
Bp=g.Bpxy
hthe=g.hthe

Dpsi=dblarr(nx,ny)
Dpsi[0,*]=psixy[1,*]-psixy[0,*]
Dpsi[nx-1,*]=psixy[nx-1,*]-psixy[nx-2,*]
FOR i=1,nx-2 DO BEGIN
  Dpsi[i,*]=(psixy[i+1,*]-psixy[i-1,*])/2
ENDFOR

Ddcp1=dblarr(nt)
Ddcp2=dblarr(nt)
Ddcp3=dblarr(nt)
Tp01=0.
Tp02=0.
Tp03=0.

FOR t=0,nt-1 DO BEGIN
 Ddcp3[t]=2.0*PI*total(mydcp[xmin:xmax,*,t]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/Bp[xmin:xmax,*])
 Ddcp2[t]=total(mydcp[xmin:xmax,*,t]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/(R[xmin:xmax,*]*Bp[xmin:xmax,*]))
 Ddcp1[t]=total(mydcp[xmin:xmax,yind,t]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind])) 
ENDFOR

  Tp03=2.0*PI*total(myp0[xmin:xmax,*]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/Bp[xmin:xmax,*])
  Tp02=total(myp0[xmin:xmax,*]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/(R[xmin:xmax,*]*Bp[xmin:xmax,*]))
  Tp01=total(myp0[xmin:xmax,yind]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind]))

s1=dblarr(nt)
s2=dblarr(nt)
s3=dblarr(nt)
E_loss=dblarr(nt)

s1=-Ddcp1/Tp01   ;1D elm size
s2=-Ddcp2/Tp02   ;2D elm size
s3=-Ddcp3/Tp03   ;3D elm size

E_loss=-Ddcp3*(0.5*Bbar*Bbar/MU0)    ;energy loss, unit J
E_total=Tp03*(0.5*Bbar*Bbar/MU0)     ;total energy, unit J

elmsize={s1:s1,s2:s2,s3:s3,E_loss:E_loss,E_total:E_total}

RETURN, elmsize
  
END
