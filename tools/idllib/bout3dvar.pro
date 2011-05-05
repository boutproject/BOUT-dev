; Converts 3D variables to and from BOUT++'s Fourier representation
; which is used for input grids.
;
; This method is used to allow the Z resolution to be varied
; without changing the input grid file.
;
; Usage:
; =====
;
; v = bout3Dvar(var)      Converts var[x,y,z] into v[x,y,f]
;                         where f are the Fourier components
; 
; By default the number of frequencies kept is nz/2
;
; var = bout3Dvar(v, /reverse)   Converts from v[x,y,f] to var[x,y,z]
; 
; By default /reverse uses nz = 2*(f-1) (so even number), but this
; can be specified manually using the nz=... keyword
;

FUNCTION bout3dvar, var, reverse=reverse, nf=nf, nz=nz
  s = SIZE(var, /dim)
  IF N_ELEMENTS(s) NE 3 THEN BEGIN
     PRINT, "Error in bout3dvar: input variable must be 3D"
     RETURN, var
  ENDIF
  
  nx = s[0]
  ny = s[1]
  
  IF NOT KEYWORD_SET(reverse) THEN BEGIN
    ; Convert from [x,y,z] to [x,y,f]
    nz = s[2]
    IF NOT KEYWORD_SET(nf) THEN nf = nz/2
     
    result = FLTARR(nx, ny, 2*nf + 1)
    
    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        fvals = FFT(var[x,y,*])
        result[x,y,0] = REAL_PART(fvals[0])
        FOR f=0,MIN([nf,nz/2])-1 DO BEGIN
          result[x,y,1+2*f] = REAL_PART(fvals[f+1])
          result[x,y,2+2*f] = IMAGINARY(fvals[f+1])
        ENDFOR
      ENDFOR
    ENDFOR
  ENDIF ELSE BEGIN
    ; Converting from frequency to real space
    nf = (s[2]-1)/2
    IF NOT KEYWORD_SET(nz) THEN nz = s[2]-1
    
    nf2 = nz/2
    fvals = CMPLXARR(nf2)
    
    result = FLTARR(nx, ny, nz)
    
    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        fvals[0] = var[x,y,0] ; DC component
        FOR f=0,nf-1 DO BEGIN
          fvals[1+f] = COMPLEX(var[x,y,1+2*f], var[x,y,2+2*f])
          fvals[nf2-1-f] = CONJUGATE(fvals[1+f])
        ENDFOR
        result[nx, ny, *] = REAL_PART(FFT(fvals,/reverse))
      ENDFOR
    ENDFOR
  ENDELSE
  
  RETURN, result
END

