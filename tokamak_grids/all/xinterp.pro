
FUNCTION xinterp, var, nx2
  s = SIZE(var)

  IF s[0] EQ 0 THEN BEGIN
      ; just a number: No interpolation
      v2 = var
  ENDIF ELSE IF s[0] EQ 1 THEN BEGIN
      s = SIZE(var, /dimensions)
      nx = s[0]
      
      newx = DINDGEN(nx2)*(DOUBLE(nx)-1.) / (DOUBLE(nx2)-1.)
      v2 = INTERPOL(var, DINDGEN(nx), newx, /QUAD)
  ENDIF ELSE BEGIN

      s = SIZE(var, /dimensions)
      nx = s[0]
      ny = s[1]
      
      newx = DINDGEN(nx2)*(DOUBLE(nx)-1.) / (DOUBLE(nx2)-1.)
      
      v2 = DBLARR(nx2, ny)
      FOR i=0, ny-1 DO BEGIN
          v2[*,i] = INTERPOL(REFORM(var[*,i]), DINDGEN(nx), newx, /QUAD)
      ENDFOR
  ENDELSE
  RETURN, v2
END
