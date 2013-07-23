FUNCTION volumeint,var2d,uedge,xmin=xmin,xmax=xmax
  
Nparameter=N_Params()
IF Nparameter NE 2 THEN BEGIN
  PRINT,"lack of parameters"
  RETURN, 0
ENDIF



var=var2d
g=uedge

PI = 3.1415926
MU0 = 4.0e-7*PI

s=size(var)

IF s[0] NE 3 THEN BEGIN
  PRINT,"var should be 3D(x,y,t)"
ENDIF

nx=s[1]
ny=s[2]
nt=s[3]

IF NOT KEYWORD_SET(xmin) THEN xmin=0
IF NOT KEYWORD_SET(xmax) THEN xmax=nx-1

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

result=dblarr(nt)

FOR t=0,nt-1 DO BEGIN
 result[t]=total(var[xmin:xmax,*,t]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/Bp[xmin:xmax,*])
ENDFOR


RETURN, result
  
END
