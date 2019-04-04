; 2D Discrete Cosine Transform 
; 
; Author: Maxim Umansky, LLNL
;
; 18-02-2010 Ben Dudson <bd512@york.ac.uk>
;      * Modified to accept non-square arrays
;      * Speed up using matrix operations instead of FOR loops
;
; 2011/06/02 Ilon Joseph
;       * Double precision for forward and inverse transformation
;       * Matrix operations for inverse
;       * Forward # Inverse = 1 with errors on the order of 1.0d-14 of maximum

pro dct2dslow, sig, fsig, inverse=inverse
;
;-calculate 2D discrete cosine transform
;----------------------------------------

IF NOT KEYWORD_SET(INVERSE) THEN BEGIN ;---direct transform---

  s = SIZE(sig, /DIMENSION)

  nx=s[0]
  ny=s[1]

  fsig=dblarr(nx,ny)
  
  for iu=0,Nx-1 do begin
    for jv=0,Ny-1 do begin
      
      fsig[iu,jv] = TOTAL( double(sig) * $
        ( COS(iu*!DPI*(2*dindgen(Nx)+1)/(2*Nx)) # COS(jv*!DPI*(2*dindgen(Ny)+1)/(2*Ny))) )
      
    endfor
  endfor

  fsig *= 2.D/SQRT(double(Nx*Ny))
  fsig[0,*] *= SQRT(0.5d0)
  fsig[*,0] *= SQRT(0.5d0)

ENDIF ELSE BEGIN ;---inverse transform---
  
  s = SIZE(fsig, /DIMENSION)

  nx=s[0]
  ny=s[1]

  dsig=double(fsig)
  dsig[0,*] *= SQRT(0.5d0)
  dsig[*,0] *= SQRT(0.5d0)

  sig=dblarr(nx,ny)

  for ix=0,Nx-1 do begin
   for jy=0,Ny-1 do begin

       sig[ix,jy]= TOTAL( dsig * $
        ( COS(dindgen(Nx)*!DPI*(2*ix+1)/(2*Nx)) # COS(dindgen(Ny)*!DPI*(2*jy+1)/(2*Ny)) ) )

   endfor
  endfor

  sig *= 2.D/SQRT(double(NX*NY))

ENDELSE
 

;
;
;
end
