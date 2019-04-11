; Integrate an array using the trapezoid rule

FUNCTION INT_TRAPEZOID, x, var
  result = 0.0D
  n = N_ELEMENTS(var)
  FOR i=1, n-1 DO BEGIN
    result = result + 0.5*(x[i] - x[i-1])*(var[i] + var[i-1])
  ENDFOR

  RETURN, result
END
