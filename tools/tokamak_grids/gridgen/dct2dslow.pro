; 2D Discrete Cosine Transform 
; 
; Author: Maxim Umansky, LLNL
;
; 18-02-2010 Ben Dudson <bd512@york.ac.uk>
;      * Modified to accept non-square arrays
;      * Speed up using matrix operations instead of FOR loops

pro dct2dslow, sig, fsig, inverse=inverse
;
;-calculate 2D discrete cosine transform
;----------------------------------------

IF NOT KEYWORD_SET(INVERSE) THEN BEGIN ;---direct transform---

  s = SIZE(sig, /DIMENSION)

  nx=s[0]
  ny=s[1]

  fsig=fltarr(nx,ny)
  
  for iu=0,Nx-1 do begin
    for jv=0,Ny-1 do begin
      
      if (iu eq 0) then cu=0.707107 else cu=1.
      if (jv eq 0) then cv=0.707107 else cv=1.
      
      sum = TOTAL( sig * ( COS(iu*!PI*(2*findgen(Nx)+1)/(2*Nx)) # COS(jv*!PI*(2*findgen(ny)+1)/(2*Ny))) )
      
      fsig[iu,jv]=SQRT(2./nx)*SQRT(2./ny)*cv*cu*sum 
      
    endfor
  endfor
  
ENDIF ELSE BEGIN ;---inverse transform---
  
  s = SIZE(fsig, /DIMENSION)

  nx=s[0]
  ny=s[1]


  sig=fltarr(nx,ny)

  for ix=0,Nx-1 do begin
   for jy=0,Ny-1 do begin

        sum=0.0

        for iu=0,Nx-1 do begin
           for jv=0,Ny-1 do begin
           ;     
           if (iu eq 0) then cu=0.707107 else cu=1.
           if (jv eq 0) then cv=0.707107 else cv=1.

            sum=sum + cv*cu*fsig[iu,jv]*cos(jv*!PI*(2*jy+1)/(2*Ny))*$
                                         cos(iu*!PI*(2*ix+1)/(2*Nx))

           ;
           endfor
        endfor

       sig[ix,jy]=SQRT(2./nx)*SQRT(2./ny)*sum ; NX only if nx=ny?
   
     ;STOP
   endfor
  endfor


ENDELSE
 

;
;
;
end
