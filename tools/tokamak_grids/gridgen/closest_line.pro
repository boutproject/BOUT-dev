; Find the closest contour line to a given point
FUNCTION closest_line, info, xy, ri, zi, mind=mind
  mind = MIN( (xy[0,info[0].offset:(info[0].offset+info[0].n-1)] - ri)^2 + $
              (xy[1,info[0].offset:(info[0].offset+info[0].n-1)] - zi)^2 )
  ind = 0
  FOR i=1, N_ELEMENTS(info)-1 DO BEGIN
    d = MIN( (xy[0,info[i].offset:(info[i].offset+info[i].n-1)] - ri)^2 + $
             (xy[1,info[i].offset:(info[i].offset+info[i].n-1)] - zi)^2 )
    IF d LT mind THEN BEGIN
      mind = d
      ind = i
    ENDIF
  ENDFOR
  RETURN, ind
END

