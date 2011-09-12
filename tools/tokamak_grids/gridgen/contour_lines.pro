; Return the contour lines of a given set of levels
; Just a small wrapper around the built-in CONTOUR
PRO contour_lines, z, x, y, levels=levels, $
                   path_info=path_info, path_xy=path_xy, $
                   _extra=_extra

  old_x = !x
  old_y = !y
  IF N_ELEMENTS(y) NE N_ELEMENTS(z) THEN BEGIN
    CONTOUR, z, levels=levels, $
      path_info=path_info, path_xy=path_xy, $
      closed=0, /PATH_DATA_COORDS, _extra=_extra
  ENDIF ELSE BEGIN
    CONTOUR, z, x, y, levels=levels, $
      path_info=path_info, path_xy=path_xy, $
      closed=0, /PATH_DATA_COORDS, _extra=_extra
  ENDELSE
  !x = old_x
  !y = old_y
END
