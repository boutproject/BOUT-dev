
FUNCTION get_float, prompt
  f = 0.0
  CATCH, var
  READ, f, prompt=prompt
  CATCH, /CANCEL
  RETURN, f
END
