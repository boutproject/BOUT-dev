
FUNCTION get_integer, prompt
  CATCH
  i = 0
  READ, i, prompt=prompt
  CATCH, /cancel
  RETURN, i
END
