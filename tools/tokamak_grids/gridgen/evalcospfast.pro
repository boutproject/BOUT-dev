function EvalCosPfast, fsig, x0=x0,y0=y0
;
;-evaluate the cos polinomial at point x0,y0
;and the partial derivatives
;
;
; Author: Maxim Umansky, LLNL
;
; 29-02-2010 Ben Dudson <bd512@york.ac.uk>
;     * Extended to non-square arrays
;--------------------------------------------

  s = SIZE(fsig, /dimens)
  nx=s[0]
  ny=s[1]

  cuvec=DBLARR(nx)+1.D
  cuvec[0]=1.D/SQRT(2.D)

  cvvec=DBLARR(ny)+1.D
  cvvec[0]=1.D/SQRT(2.D)
  
  uvec=COS(!DPI*findgen(nx)*(x0+0.5D)/nx)
  uvex=(-findgen(nx)*!DPI/nx)*SIN(!DPI*findgen(nx)*(x0+0.5D)/nx)
  
  vvec=COS(!DPI*findgen(ny)*(y0+0.5D)/ny)
  vvey=(-findgen(ny)*!DPI/ny)*SIN(!DPI*findgen(ny)*(y0+0.5D)/ny)
  
  ;-value
  res=SQRT(2.D/nx)*SQRT(2.D/ny) * TOTAL(((cuvec # cvvec) * fsig) * (uvec # vvec))  
  
  ;d/dx
  rex=SQRT(2.D/nx)*SQRT(2.D/ny) * TOTAL(((cuvec # cvvec) * fsig) * (uvex # vvec))  
  
  ;d/dy
  rey=SQRT(2.D/nx)*SQRT(2.D/ny) * TOTAL(((cuvec # cvvec) * fsig) * (uvec # vvey))  
  
;
;
;
return, [res,rex,rey]
end

