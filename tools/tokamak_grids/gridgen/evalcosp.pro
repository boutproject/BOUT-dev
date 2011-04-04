function EvalCosP, fsig, x0=x0,y0=y0
;
; evaluate the cos polinomial at point x0,y0
; and the partial derivatives
;
; Author: Maxim Umansky, LLNL
;
; 20-01-2010 Ben Dudson <bd512@york.ac.uk>
;      * Extended to calculate second derivatives
;
; 29-02-2010 Ben Dudson <bd512@york.ac.uk>
;      * Extended to non-square arrays
;--------------------------------------------
  
  s = SIZE(fsig, /dimens)
  nx=s[0]
  ny=s[1]

  sum=0.0
  ; First derivatives
  sumx=0.0
  sumy=0.0
  ; Second derivatives
  sumxx = 0.0
  sumyy = 0.0
  sumxy = 0.0
  
  for iu=0,Nx-1 do begin
    for jv=0,Ny-1 do begin
      ;     
      if (iu eq 0) then cu=0.707107 else cu=1.
      if (jv eq 0) then cv=0.707107 else cv=1.
      
      sum=sum + cv*cu*fsig[iu,jv]*$
        COS(jv*!PI*(2*y0+1)/(2*Ny))*COS(iu*!PI*(2*x0+1)/(2*Nx))
      
      sumx=sumx + cv*cu*fsig[iu,jv]*$
        COS(jv*!PI*(2*y0+1)/(2*Ny))*SIN(iu*!PI*(2*x0+1)/(2*Nx))*$
        (-iu*!PI/Nx)
      
      sumy=sumy + cv*cu*fsig[iu,jv]*$
        SIN(jv*!PI*(2*y0+1)/(2*Ny))*COS(iu*!PI*(2*x0+1)/(2*Nx))*$
        (-jv*!PI/Ny)
      
      sumxx = sumxx - cv*cu*fsig[iu,jv]*$
        COS(jv*!PI*(2*y0+1)/(2*Ny))*COS(iu*!PI*(2*x0+1)/(2*Nx))*$
        (iu*!PI/Nx)^2
      
      sumyy = sumyy - cv*cu*fsig[iu,jv]*$
        COS(jv*!PI*(2*y0+1)/(2*Ny))*COS(iu*!PI*(2*x0+1)/(2*Nx))*$
        (jv*!PI/Ny)^2
      
      sumxy = sumxy + cv*cu*fsig[iu,jv]*$
        SIN(jv*!PI*(2*y0+1)/(2*Ny))*SIN(iu*!PI*(2*x0+1)/(2*Nx))*$
        (iu*!PI/Nx)*(jv*!PI/Ny)
      ;
    endfor
  endfor
  
  res=SQRT(2./Nx)*SQRT(2./Ny)*[sum, sumx,sumy, sumxx,sumyy,sumxy]   
;
;
;
return, res
end
