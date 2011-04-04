
FUNCTION zshift, var, shift, period=period

  IF NOT KEYWORD_SET(period) THEN period = 1.0

  var = REFORM(var)
  shift = REFORM(shift)

  s = SIZE(var)
  s2 = SIZE(shift)

  newvar = 0

  IF (s[0] EQ 3) AND (s2[0] EQ 1) THEN BEGIN
      PRINT, "Shifting an XZT variable"
      
      nx = s[1]
      nz = s[2]
      nt = s[3]
      
      dz = 2.0*!PI / (period * FLOAT(nz-1))
      
      newvar = FLTARR(nx, nz, nt)

      FOR x=0, nx-1 DO BEGIN
          offset = shift[x] / dz
          
          FOR z=0, nz-1 DO BEGIN
              zpos = (((z + offset) MOD (nz-2)) + (nz-2)) MOD (nz-2)
              
              FOR t=0, nt-1 DO BEGIN
                  newvar[x,z,t] = INTERPOL(REFORM(var[x,*,t]), FINDGEN(nz), zpos, /SPLINE)
              ENDFOR
          ENDFOR
      ENDFOR
      
  ENDIF ELSE IF (s[0] EQ 2) AND (s2[0] EQ 1) THEN BEGIN
     PRINT, "Shifting an XZ variable"

     nx = s[1]
     nz = s[2]
     dz = 2.0*!PI / (period * FLOAT(nz-1))
      
     newvar = FLTARR(nx, nz)
     
     FOR x=0, nx-1 DO BEGIN
        offset = shift[x] / dz
        
        FOR z=0, nz-1 DO BEGIN
           zpos = (((z + offset) MOD (nz-2)) + (nz-2)) MOD (nz-2)
           
           newvar[x,z] = INTERPOL(REFORM(var[x,*]), FINDGEN(nz), zpos, /SPLINE)
        ENDFOR
     ENDFOR
  ENDIF ELSE BEGIN
      PRINT, "Don't know what to do"
  ENDELSE
  
  RETURN, newvar
END
