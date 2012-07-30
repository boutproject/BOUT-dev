; Tests if a given integer is prime

FUNCTION is_prime, val
  n = LONG(val)

  IF n LE 3 THEN RETURN, 1

  i = 2L
  sn = SQRT(n)
  
  REPEAT BEGIN
    IF n MOD i EQ 0 THEN RETURN, 0
    i = i + 1L
  ENDREP UNTIL i GT sn
  RETURN, 1
END

