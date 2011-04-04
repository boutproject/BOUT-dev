; Extract all matches to pattern contained in str
; 
; Designed for extracting numbers from FORTRAN-formatted
; text output files.
;
; Example: To extract all numbers from a string, use
;
;   line = '-0.341859439E+01-0.341601279E+01-0.341354986E+01'
;   pattern = '[+-]?([0-9])+($| |(\.[0-9]+([eE][+-]?[0-9]*)?))'
;   data = regex_extract(line, pattern)
; 
;   data is now an array of strings
;      ['-0.341859439E+01', '-0.341601279E+01', '-0.341354986E+01']
;
; Ben Dudson, University of York, Feb 2010

FUNCTION regex_extract, line, pattern, nmatch=nmatch
  result = 0
  nmatch = 0
  str = line
  REPEAT BEGIN
    ind = STREGEX(str, pattern)
    IF ind GE 0 THEN BEGIN
      s = STREGEX(str, pattern, /extract)
      IF ind + STRLEN(s) LT STRLEN(str) THEN BEGIN
        str = STRMID(str, ind+STRLEN(s))
      ENDIF ELSE ind = -1
      IF nmatch EQ 0 THEN result = [s] ELSE result = [result, s]
      nmatch = nmatch + 1
    ENDIF
  ENDREP UNTIL ind LT 0
  RETURN, result
END

