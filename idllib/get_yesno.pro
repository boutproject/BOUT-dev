; asks the user for a yes/no response

FUNCTION get_yesno, prompt, gui=gui, _extra=_extra
  IF KEYWORD_SET(gui) THEN BEGIN
    r = DIALOG_MESSAGE(prompt, /question, _extra=_extra)
    RETURN, STRCMP(r, "Yes", /fold_case)
  ENDIF ELSE BEGIN
    s = "yes"
    REPEAT BEGIN
      READ, s, prompt=prompt
      s = STRLOWCASE(s)
      IF (s EQ "y") OR (s EQ "yes") THEN BEGIN
        RETURN,1
      ENDIF ELSE IF (s EQ "n") OR (s EQ "no") THEN BEGIN
        RETURN, 0
      ENDIF 
    ENDREP UNTIL 0
  ENDELSE
END
