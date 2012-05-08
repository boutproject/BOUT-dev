;
; Read .equ tokamak R-Z equilibrium file
; 
; 
; Ben Dudson, University of York, Feb 2012
;

FUNCTION toFloat, text
  num = 0.0
  READS, text, num
  RETURN, num
END

FUNCTION read_equ, file
  OPENR, fid, file, /GET_LUN, error=errid
  IF errid THEN BEGIN
    PRINT, "Error whilst reading '"+file+"' :" + !ERROR_STATE.MSG
    RETURN, 0
  ENDIF

  btf = 0.0
  rf = 0.0
  
  REPEAT BEGIN
     ; read the next line
     line=' '
     READF, fid, line

     ; strip spaces
     line = STRTRIM(line, 2)
     IF STRLEN(line) EQ 0 THEN CONTINUE
     
     CASE line OF
       'r(1:jm);': r = read_1d(jm, fid=fid)
       'z(1:km);': z = read_1d(km, fid=fid)
       '((psi(j,k)-psib,j=1,jm),k=1,km)': BEGIN
           psi = read_2d(jm, km, fid=fid)
        END
        ELSE: BEGIN
          IF STREGEX(line, ':=') GE 0 THEN BEGIN
            ; contains ':='
            s = STRSPLIT(line, ':=', /extract)

            left  = STRTRIM(s[0],2)
            right = STRTRIM(s[1],2)

          ENDIF ELSE IF STREGEX(line, '=') GE 0 THEN BEGIN
            ; contains '='
            s = STRSPLIT(line, '=', /extract)

            left  = STRLOWCASE(STRTRIM(s[0],2))
            right = STRTRIM(s[1],2)

            PRINT, "Equality: '"+left+"' = '"+right+"'"
            CASE left OF
               'jm': jm = FIX(right)
               'km': km = FIX(right)
               'psib': psib = toFloat(right)
               'btf': btf = toFloat(right)
               'rf': rf = toFloat(right)
                ELSE: BEGIN
                  PRINT, "Unrecognised variable: "+left
                END
            ENDCASE
          ENDIF ELSE BEGIN
            PRINT, "Unrecognised line: '"+ line+"'"
            STOP
          ENDELSE
        END
      ENDCASE
     
  ENDREP UNTIL EOF(fid)

  FREE_LUN, fid
  
  result = {nx:jm, ny:km, $
            r:r, z:z, $
            psi:psi, $
            fpol:(btf * rf)}
  
  RETURN, result
END
