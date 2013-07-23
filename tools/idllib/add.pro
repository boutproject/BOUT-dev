FUNCTION add, var1,var2,complex=complex

s1=size(var1)
s2=size(var2)

IF s1[0] EQ 3 THEN BEGIN
t1=s1[3]
t2=s2[3]
nx=s1[1]
ny=s1[2]

IF KEYWORD_SET(complex) THEN result=complexarr(nx,ny,t1+t2) ELSE result=dblarr(nx,ny,t1+t2)
for t=0, t1-1 do result[*,*,t]=var1[*,*,t]
for t=0, t2-1 do result[*,*,t1+t]=var2[*,*,t]

ENDIF ELSE IF s1[0] EQ 4 THEN BEGIN

t1=s1[4]
t2=s2[4]
nx=s1[1]
ny=s1[2]
nz=s1[3]

IF KEYWORD_SET(complex) THEN result=complexarr(nx,ny,nz,t1+t2) ELSE result=dblarr(nx,ny,nz,t1+t2)
for t=0, t1-1 do result[*,*,*,t]=var1[*,*,*,t]
for t=0, t2-1 do result[*,*,*,t1+t]=var2[*,*,*,t]
 
ENDIF ELSE IF s1[0] EQ 2 THEN BEGIN
t1=s1[2]
t2=s2[2]
nx=s1[1]

IF KEYWORD_SET(complex) THEN result=complexarr(nx,t1+t2) ELSE result=dblarr(nx,t1+t2)
for t=0, t1-1 do result[*,t]=var1[*,t]
for t=0, t2-1 do result[*,t1+t]=var2[*,t]
ENDIF

return, result

end
