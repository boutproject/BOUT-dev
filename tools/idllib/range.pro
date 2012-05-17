; A useful little function to generate ranges

FUNCTION range, first, last
  IF first LT last THEN BEGIN
    RETURN, first + INDGEN(last - first + 1)
  ENDIF ELSE BEGIN
    RETURN, last + REVERSE(INDGEN(first - last + 1))
  ENDELSE
END

