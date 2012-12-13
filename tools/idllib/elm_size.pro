FUNCTION elm_size,dcp,p0,uedge,xmin=xmin,xmax=xmax,dim=dim,yind=yind

Nparameter=N_Params()
IF Nparameter NE 3 THEN BEGIN
  PRINT,"lack of parameters"
  RETURN, 0
ENDIF

IF NOT KEYWORD_SET(xmin) THEN xmin=0
IF NOT KEYWORD_SET(xmax) THEN xmax=327  
IF NOT KEYWORD_SET(dim)  THEN dim=1
IF NOT KEYWORK_SET(yind) THEN yind=32  ; outer midplane location

mydcp=dcp
myp0=p0
g=uedge

s=size(mydcp)

IF s[0] NE 3 THEN BEGIN
  PRINT,"dcp should be 3D(x,y,t)"
ENDIF

nx=s[1]
ny=s[2]
nt=s[3]

theta=g.dy     ;using correct poloidal angle
psixy=g.psixy
R=g.rxy
Bp=g.Bpxy
hthe=g.hthe

Dpsi=dblarr(nx,ny)
Dtheta=dblarr(nx,ny)

Dpsi[0,*]=(psixy[1,*]-psixy[0,*])/2
Dpsi[nx-1,*]=(psixy[nx-1,*]-psixy[nx-2,*])/2

FOR i=1,nx-2 DO BEGIN
  Dpsi[i,*]=(psixy[i+1,*]-psixy[i-1,*])/2
ENDFOR

Dtheta[*,0]=(2*3.1415926+theta[*,1]-theta[*,ny-1])/2
Dtheta[*,ny-1]=(2*3.1415926-theta[*,ny-2])/2

FOR j=1,ny-2 DO BEGIN
  Dtheta[*,j]=(theta[*,j+1]-theta[*,j-1])/2
ENDFOR

Ddcp=dblarr(nt)
Tp0=0.

FOR t=0,nt-1 DO BEGIN
IF dim EQ 3 THEN BEGIN           ;3D integral gives thermal energy unit
 Ddcp[t]=total(mydcp[xmin:xmax,*,t]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/Bp[xmin:xmax,*])
ENDIF

IF dim EQ 2 THEN BEGIN          ;2D integral gives the relative energy loss on a cross-section
 Ddcp[t]=total(mydcp[xmin:xmax,*,t]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/(R[xmin:xmax,*]*Bp[xmin:xmax,*]))
ENDIF

IF dim EQ 1 THEN BEGIN          ;1D integral gives the relateive energy loss at outer midplane
 Ddcp[t]=total(mydcp[xmin:xmax,yind,t]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind])) 
ENDIF
ENDFOR

IF dim EQ 3 THEN BEGIN
  Tp0=total(myp0[xmin:xmax,*]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/Bp[xmin:xmax,*])
ENDIF
IF dim EQ 2 THEN BEGIN
  Tp0=total(myp0[xmin:xmax,*]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/(R[xmin:xmax,*]*Bp[xmin:xmax,*]))
ENDIF
IF dim EQ 1 THEN BEGIN
  Tp0=total(myp0[xmin:xmax,yind]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind]))
ENDIF

elmsize=-Ddcp/Tp0

RETURN, elmsize
  
END
