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

  cuvec=DBLARR(nx)+1.
  cuvec[0]=1./SQRT(2.)

  cvvec=DBLARR(ny)+1.
  cvvec[0]=1./SQRT(2.)
  
  uvec=COS(!PI*findgen(nx)*(x0+0.5)/nx)
  uvex=(-findgen(nx)*!PI/nx)*SIN(!PI*findgen(nx)*(x0+0.5)/nx)
  
  vvec=COS(!PI*findgen(ny)*(y0+0.5)/ny)
  vvey=(-findgen(ny)*!PI/ny)*SIN(!PI*findgen(ny)*(y0+0.5)/ny)
  
  ;-value
  res=SQRT(2./nx)*SQRT(2./ny) * TOTAL(((cuvec # cvvec) * fsig) * (uvec # vvec))  
  
  ;d/dx
  rex=SQRT(2./nx)*SQRT(2./ny) * TOTAL(((cuvec # cvvec) * fsig) * (uvex # vvec))  
  
  ;d/dy
  rey=SQRT(2./nx)*SQRT(2./ny) * TOTAL(((cuvec # cvvec) * fsig) * (uvec # vvey))  
  
;
;
;
return, [res,rex,rey]
end

