FUNCTION elm_size_yind,dcp,p0,uedge,geom=geom,xmin=xmin,xmax=xmax, ymin=ymin, ymax=ymax, yind=yind, dim=dim

Nparameter=N_Params()
IF Nparameter NE 3 THEN BEGIN
  PRINT,"lack of parameters"
  RETURN, 0
ENDIF

   PRINT,1,xmin, xmax,ymin,ymax, yind
IF NOT KEYWORD_SET(xmin) THEN xmin=0
IF NOT KEYWORD_SET(xmax) THEN xmax=327  
;IF NOT KEYWORD_SET(yind) THEN yind=75  
IF NOT KEYWORD_SET(dim)  THEN dim=1
IF NOT KEYWORD_SET(geom)  THEN geom=circle

PRINT,2,xmin, xmax,ymin,ymax, yind

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

theta=g.pol_angle
psixy=g.psixy
R=g.rxy
Bp=g.Bpxy
hthe=g.hthe

Dpsi=dblarr(nx,ny)
Dtheta=dblarr(nx,ny)

Dpsi[0,*]=(psixy[1,*]-psixy[0,*])
Dpsi[nx-1,*]=(psixy[nx-1,*]-psixy[nx-2,*])

FOR i=1,nx-2 DO BEGIN
  Dpsi[i,*]=(psixy[i+1,*]-psixy[i-1,*])/2
ENDFOR

Ddcp=dblarr(nt)
Tp0=0.

; for circlular geometry
IF geom eq "circle" THEN BEGIN
   Dtheta[*,0]=(2*3.1415926+theta[*,1]-theta[*,ny-1])/2
   Dtheta[*,ny-1]=(2*3.1415926-theta[*,ny-2])/2

   FOR j=1,ny-2 DO BEGIN
      Dtheta[*,j]=(theta[*,j+1]-theta[*,j-1])/2
   ENDFOR

   FOR t=0,nt-1 DO BEGIN
      IF dim EQ 2 THEN BEGIN
         Ddcp[t]=total(mydcp[xmin:xmax,*,t]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/(R[xmin:xmax,*]*Bp[xmin:xmax,*]))
      ENDIF ELSE BEGIN
         Ddcp[t]=total(mydcp[xmin:xmax,yind,t]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind])) 
      ENDELSE
   ENDFOR

   IF dim EQ 2 THEN BEGIN
      Tp0=total(myp0[xmin:xmax,*]*hthe[xmin:xmax,*]*Dtheta[xmin:xmax,*]*Dpsi[xmin:xmax,*]/(R[xmin:xmax,*]*Bp[xmin:xmax,*]))
   ENDIF ELSE BEGIN
      Tp0=total(myp0[xmin:xmax,yind]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind]))
   ENDELSE
ENDIF

   PRINT,3,xmin, xmax,ymin,ymax, yind

; for divertor geometry
IF geom eq "snl" THEN BEGIN
   PRINT,"geom=", geom
   xmax1 = xmax
   IF (yind eq 0) OR (yind eq ny-1) THEN BEGIN
      xmin1=0 
      xmax1=nx-1

      FOR t=0,nt-1 DO BEGIN
         Ddcp[t]=-total(mydcp[xmin1:xmax1,yind,t]*Dpsi[xmin1:xmax1,yind]/(R[xmin1:xmax1,yind]*Bp[xmin1:xmax1,yind])) 
      ENDFOR
   ENDIF
   
   Dtheta[*,0]=(2*3.1415926+theta[*,ymin]-theta[*,ymax])
   Dtheta[*,ymax]=(2*3.1415926-theta[*,ymax])

   FOR j=ymin+1,ymax-1 DO BEGIN
      Dtheta[*,j]=(theta[*,j+1]-theta[*,j-1])/2
   ENDFOR

   IF yind eq 75 THEN BEGIN
      FOR t=0,nt-1 DO BEGIN
         IF dim EQ 2 THEN BEGIN
            Ddcp[t]=total(mydcp[xmin:xmax,ymin:ymax,t]*hthe[xmin:xmax,ymin:ymax]*Dtheta[xmin:xmax,ymin:ymax]*Dpsi[xmin:xmax,ymin:ymax]/(R[xmin:xmax,ymin:ymax]*Bp[xmin:xmax,ymin:ymax]))
         ENDIF ELSE IF dim EQ 1 THEN BEGIN
            Ddcp[t]=total(mydcp[xmin:xmax,yind,t]*Dpsi[xmin:xmax,yind]/(R[xmin:xmax,yind]*Bp[xmin:xmax,yind])) 
         ENDIF
      ENDFOR
   ENDIF
   PRINT,xmin, xmax,ymin,ymax,yind

   Tp0=total(myp0[xmin:xmax,ymin:ymax]*hthe[xmin:xmax,ymin:ymax]*Dtheta[xmin:xmax,ymin:ymax]*Dpsi[xmin:xmax,ymin:ymax]/(R[xmin:xmax,ymin:ymax]*Bp[xmin:xmax,ymin:ymax]))
   PRINT,"Tp0=", Tp0

ENDIF

elmsize=-Ddcp/Tp0

RETURN, elmsize
  
END
