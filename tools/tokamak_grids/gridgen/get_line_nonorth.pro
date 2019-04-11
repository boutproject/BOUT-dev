FUNCTION get_line_nonorth, interp_data, R, Z, ri0, zi0, fto, npt=npt, vec_up=vec_up, weight_up=weight_up, vec_down=vec_down, weight_down=weight_down
  IF NOT KEYWORD_SET(npt) THEN npt=10
  ; Get starting f
  
  local_gradient, interp_data, ri0, zi0, status=status, f=ffrom
  IF status THEN BEGIN
    PRINT, "WARNING: get_line starting location "+STR([ri0,zi0])+" is out of domain"
    RETURN, [[ri0,ri0],[zi0,zi0]]
  ENDIF

  IF NOT KEYWORD_SET(weight_up) THEN weight_up = 0.D
  IF NOT KEYWORD_SET(weight_down) THEN weight_down = 0.D
  rixpt = DBLARR(npt+1)
  zixpt = rixpt
  rixpt[0] = ri0
  zixpt[0] = zi0
  FOR j=0, npt-1 DO BEGIN
    d = DOUBLE(j+1)/DOUBLE(npt)
    ftarg = d*fto + (1.0D - d)*ffrom
    follow_gradient_nonorth, interp_data, R, Z, rixpt[j], zixpt[j], $
      ftarg, rinext, zinext, vec_up=vec_up, weight_up=weight_up, $
      vec_down=vec_down, weight_down=weight_down
    rixpt[j+1] = rinext
    zixpt[j+1] = zinext
  ENDFOR
  
  RETURN, [[rixpt], [zixpt]]
END
