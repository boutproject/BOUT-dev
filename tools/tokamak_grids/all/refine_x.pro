; refine a radial grid 
;
; Can do the following:
; - Make sure grid-points are located at rational surfaces
; - put more grid-points near rational surfaces
; - Refine grid so that the toroidal shift between surfaces
;   over a poloidal step is less than one grid-point

; INPUTS
; nmax, is maximum n number for rational surfaces
; nz is the number of points in 2pi (dz = 2pi / nz)
;

FUNCTION refine_x, Rxy, Zxy, zshift, q, nmax=nmax, $
                   nz=nz, shear=shear
  ; check dimensions

  s = SIZE(Rxy, /dim)
  s2 = SIZE(Zxy, /dim)
  s3 = SIZE(zshift, /dim)
  IF N_ELEMENTS(s) NE 2 THEN BEGIN
     PRINT, "Error: Rxy must be a 2D array"
     RETURN, 0
  ENDIF

  IF MAX([s - s2, s - s3]) NE 0 THEN BEGIN
     PRINT, "Error: Rxy, Zxy and zshift must have the same dimensions"
     RETURN, 0
  ENDIF

  nx = s[0]
  ny = s[1]

  IF N_ELEMENTS(q) NE nx THEN BEGIN
     PRINT, "Error: q must have nx points"
     RETURN, 0
  ENDIF

  ; set defaults

  IF NOT KEYWORD_SET(nz) THEN nz = 128
  IF NOT KEYWORD_SET(nmax) THEN nmax = nz

  xarr = FINDGEN(nx)
  dz = 2.0*!PI / FLOAT(nz)

  ; find indices of rational surfaces

  found = 0
  FOR n=1, nmax DO BEGIN
     FOR m=1, n DO BEGIN
        nm = FLOAT(n) / FLOAT(m)
        IF (nm LE MAX(q)) AND (nm GE MIN(q)) THEN BEGIN 
           i = INTERPOL(xarr, q, nm)
           
           IF found EQ 0 THEN BEGIN
              rational = [i]
              found = 1
           ENDIF ELSE rational = [rational, i]
        ENDIF
     ENDFOR
  ENDFOR

  w = SORT(rational) ; put in ascending order
  rational = rational[w]

  ; remove duplicates (more cunning way?)


  REPEAT BEGIN
      rem = WHERE((rational[1:*] - rational[0:N_ELEMENTS(rational)-2]) LT 1.0e-4, count)
      ;STOP
      IF count GT 0 THEN BEGIN
          PRINT, "removing ", count
          ; some duplicates
          found = 0
          
          i = 0L
          REPEAT BEGIN    
              w = WHERE(rem EQ i, c)
              IF c EQ 0 THEN BEGIN
                  ; not in the remove list
                  IF found EQ 0 THEN BEGIN
                      keep = [i]
                      found = 1
                  ENDIF ELSE BEGIN
                      keep = [keep, i]
                  ENDELSE
              ENDIF
              i = i + 1L
          ENDREP UNTIL i EQ N_ELEMENTS(rational)
          

          rational = rational[keep]
      ENDIF
  ENDREP UNTIL count EQ 0
  
  IF found EQ 0 THEN BEGIN
     PRINT, "No rational surfaces found"
     STOP
  ENDIF ELSE BEGIN
     PRINT, "Number of rational surfaces: ", N_ELEMENTS(rational)
  ENDELSE
  
  ; plot location of rational surfaces

  safe_colors, /first

  PLOT, q, color=1
  FOR i=0, N_ELEMENTS(rational)-1 DO BEGIN
     OPLOT, [rational[i], rational[i]], [0, 50.], color=1
  ENDFOR

  cursor, x, y, /down

  plot, rational, color=1
  cursor, x, y, /down

  ; add a point at beginning and end
  points = [0.0, rational, nx-1]

  IF KEYWORD_SET(shear) THEN BEGIN
     ; insert grid-points where needed to ensure that toroidal
     ; shift is less than one grid-point

     pos = 0
     REPEAT BEGIN
        ; check maximum toroidal shift between two surfaces
        
        FOR y=0, ny-2 DO BEGIN
           ; toroidal shift on this flux-surface
           dz0 = INTERPOL(zshift[*,y+1], xarr, points[pos]) - $
                 INTERPOL(zshift[*,y], xarr, points[pos])

           ; toroidal shift on next surface
           dz1 = INTERPOL(zshift[*,y+1], xarr, points[pos+1]) - $
                 INTERPOL(zshift[*,y], xarr, points[pos+1])
           
           IF ABS(dz1 - dz0) / dz GT 1.0 THEN BEGIN
              ; insert a grid-point between these points

              points = [points[0:pos], $
                        (points[pos] + points[pos+1])/2.0, $
                        points[pos+1:*]]
           ENDIF
        ENDFOR
        pos = pos + 1
        
     ENDREP UNTIL pos EQ N_ELEMENTS(points) - 1
  ENDIF

  PLOT, points

  cursor, x, y, /down

  PLOT, interpol(rxy[*, FIX(ny/2)], findgen(ny), points) - rxy[0,FIX(ny/2)]

  STOP
  RETURN, points
END
