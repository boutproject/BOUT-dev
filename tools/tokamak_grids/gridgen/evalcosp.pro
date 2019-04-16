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

  sum=0.0D
  ; First derivatives
  sumx=0.0D
  sumy=0.0D
  ; Second derivatives
  sumxx = 0.0D
  sumyy = 0.0D
  sumxy = 0.0D
  
  for iu=0,Nx-1 do begin
    for jv=0,Ny-1 do begin
      ;     
      if (iu eq 0) then cu=0.707107D else cu=1.D
      if (jv eq 0) then cv=0.707107D else cv=1.D
      
      sum=sum + cv*cu*fsig[iu,jv]*$
        COS(jv*!DPI*(2*y0+1)/(2*Ny))*COS(iu*!DPI*(2*x0+1)/(2*Nx))
      
      sumx=sumx + cv*cu*fsig[iu,jv]*$
        COS(jv*!DPI*(2*y0+1)/(2*Ny))*SIN(iu*!DPI*(2*x0+1)/(2*Nx))*$
        (-iu*!DPI/Nx)
      
      sumy=sumy + cv*cu*fsig[iu,jv]*$
        SIN(jv*!DPI*(2*y0+1)/(2*Ny))*COS(iu*!DPI*(2*x0+1)/(2*Nx))*$
        (-jv*!DPI/Ny)
      
      sumxx = sumxx - cv*cu*fsig[iu,jv]*$
        COS(jv*!DPI*(2*y0+1)/(2*Ny))*COS(iu*!DPI*(2*x0+1)/(2*Nx))*$
        (iu*!DPI/Nx)^2
      
      sumyy = sumyy - cv*cu*fsig[iu,jv]*$
        COS(jv*!DPI*(2*y0+1)/(2*Ny))*COS(iu*!DPI*(2*x0+1)/(2*Nx))*$
        (jv*!DPI/Ny)^2
      
      sumxy = sumxy + cv*cu*fsig[iu,jv]*$
        SIN(jv*!DPI*(2*y0+1)/(2*Ny))*SIN(iu*!DPI*(2*x0+1)/(2*Nx))*$
        (iu*!DPI/Nx)*(jv*!DPI/Ny)
      ;
    endfor
  endfor
  
  res=SQRT(2.D/Nx)*SQRT(2.D/Ny)*[sum, sumx,sumy, sumxx,sumyy,sumxy]   
;
;
;
return, res
end
