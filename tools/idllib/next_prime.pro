; Find the next prime above given number

FUNCTION next_prime, val
  n = LONG(val)

  REPEAT BEGIN
    n = n + 1L
  ENDREP UNTIL is_prime(n)
  RETURN, n
END
