; Takes a 3D variable, and returns a 2D slice at fixed toroidal angle
; uedge is a structure containing UEDGE.GRD.PDB data
; N sets the number of times the data must be repeated for a full
; torus, e.g. n=2 is half a torus
; zangle gives the (real) toroidal angle of the result

FUNCTION get_pol_slice, var3d, uedge, n=n, zangle=zangle
  s = SIZE(var3d, /dimensions)
  IF N_ELEMENTS(s) NE 3 THEN BEGIN
      PRINT, "Error: Variable must be 3 dimensional"
      RETURN, 0
  ENDIF
  
  IF NOT KEYWORD_SET(n) THEN n = 1  ; default = full torus
  IF NOT KEYWORD_SET(zangle) THEN zangle = 0.0

  n = FIX(n) ; make sure it's an integer
  IF n LT 1 THEN n = 1
  zangle = FLOAT(zangle)

  nx = s[0]
  ny = s[1]
  nz = s[2]

  dz = 2.0*!PI / FLOAT(n*(nz-1))

  var2d = FLTARR(nx, ny)
  

  ; GET THE TOROIDAL SHIFT
  tn = TAG_NAMES(uedge)
  tn = STRUPCASE(tn)
  w = WHERE(tn EQ "QINTY", count)
  IF count GT 0 THEN BEGIN
      PRINT, "Using qinty as toroidal shift angle"
      zShift = uedge.qinty
  ENDIF ELSE BEGIN
      w = WHERE(tn EQ "ZSHIFT", count)
      IF count GT 0 THEN BEGIN
          PRINT, "Using zShift as toroidal shift angle"
          zShift = uedge.zShift
      ENDIF ELSE BEGIN
          PRINT, "ERROR: Can't find qinty or zShift variable"
          RETURN, 0
      ENDELSE
  ENDELSE


  FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
          zoffset = zShift[x, y]

          zind = (zangle - zoffset)/dz ; get the index for this angle
          ;print, x, y, zind

          z0 = ROUND(zind)

          p = zind - FLOAT(z0) ; between -0.5 and 0.5

          IF p LT 0.0 THEN BEGIN
              z0 = z0 - 1
              p = p + 1.0
          ENDIF
          
          z0 = ((z0 MOD (nz-1)) + (nz-1)) MOD (nz-1) ; get between 0 and nz-2
          
          ;PRINT, "x,y, z=", x,y, zind, z0
          
          ; for now 3-point interpolation
          
          zp = (z0 + 1) MOD (nz - 1)
          zm = (z0 - 1 + (nz-1)) MOD (nz - 1)

          var2d[x,y] = 0.5*p*(p-1.0)*var3d[x,y,zm] $
            + (1.0 - p*p)*var3d[x,y,z0] $
            + 0.5*p*(p+1.0)*var3d[x,y,zp]
          
      ENDFOR
  ENDFOR
  
  RETURN, var2d
END
