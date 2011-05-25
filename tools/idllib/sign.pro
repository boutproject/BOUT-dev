FUNCTION sign, var
  IF var GT 0. THEN RETURN, 1.
  RETURN, -1.
END
