FUNCTION total_p, dcp,p0

N_parameters=N_params()
IF N_parameters NE 2 THEN BEGIN
 PRING,"Need dcp, p0 and grid"
 RETURN,0
ENDIF

mydcp=dcp
myp0=p0
s=size(mydcp)



IF s[0] NE 3 THEN BEGIN 
  PRINT, "Error: variable must be 3D(x,y,t)"
  RETURN,0
ENDIF

x=s[1]
y=s[2]
t=s[3]

a=fltarr(x,y,t)
result=fltarr(x,y,t)

FOR i=0,t-1 DO BEGIN
     a[*,*,i]=myp0[*,*]
ENDFOR

result=a+mydcp

RETURN,result

end
