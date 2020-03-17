PRO oplot_contour, info, xy, R, Z, periodic=periodic, _extra=_extra
  ri = REFORM(xy[0,info.offset:(info.offset+info.n-1)])
  zi = REFORM(xy[1,info.offset:(info.offset+info.n-1)])
  
  IF KEYWORD_SET(periodic) THEN BEGIN
    ri = [ri, ri[0]]
    zi = [zi, zi[0]]
  ENDIF
  OPLOT, INTERPOLATE(R, ri, /DOUBLE), INTERPOLATE(Z, zi, /DOUBLE), _extra=_extra
END

