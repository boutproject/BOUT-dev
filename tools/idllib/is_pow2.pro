;; Returns 1 (true) if the given number is a power of 2

FUNCTION is_pow2, num
  ON_ERROR, 2
  
  x = LONG(num)

  RETURN, (x GT 0) AND (( (x - 1) AND x ) EQ 0)
END
