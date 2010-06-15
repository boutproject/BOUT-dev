;; Remove all except a given Z mode from the data

FUNCTION zfilter, var, n
  ON_ERROR, 2

  s = SIZE(var, /dim)

  IF n LT 1 THEN n = 1

  IF N_ELEMENTS(s) EQ 3 THEN BEGIN
    PRINT, "Filtering XYZ variable"
    nx = s[0]
    ny = s[1]
    nz = s[2]
    IF nz MOD 2 EQ 1 THEN nz = nz - 1

    IF n GT (nz/2) THEN n = nz/2
    
    var2 = FLTARR(nx, ny, nz)

    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        f = FFT(var[x,y,0:(nz-1)])
        fp = f[n]
        fm = f[nz-n]
        f = f*0.0
        f[n] = fp
        f[nz-n] = fm
        var2[x,y,*] = REAL_PART(FFT(f, /INV))
      ENDFOR
    ENDFOR
  ENDIF ELSE IF N_ELEMENTS(s) EQ 4 THEN BEGIN
    PRINT, "Filtering XYZT variable"
    
    PRINT, "Filtering XYZ variable"
    nx = s[0]
    ny = s[1]
    nz = s[2]
    IF nz MOD 2 EQ 1 THEN nz = nz - 1
    nt = s[3]

    IF n GT (nz/2) THEN n = nz/2

    var2 = FLTARR(nx, ny, nz, nt)

    FOR x=0, nx-1 DO BEGIN
      FOR y=0, ny-1 DO BEGIN
        FOR t=0, nt-1 DO BEGIN
          f = FFT(var[x,y,0:(nz-1)])
          fp = f[n]
          fm = f[nz-n]
          f = f*0.0
          f[n] = fp
          f[nz-n] = fm
          var2[x,y,*,t] = REAL_PART(FFT(f, /INV))
        ENDFOR
      ENDFOR
    ENDFOR
  ENDIF ELSE BEGIN
    PRINT, "Error: Variable must be 3 or 4-dimensional"
    RETURN, 0
  ENDELSE
  
  
  RETURN, var2
END
