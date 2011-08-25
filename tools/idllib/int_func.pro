; integrate a function, always using the maximum
; number of grid-points possible for highest accuracy
;
; Changelog
; ---------
;
; 2010-05-24 Ben Dudson <bd512@york.ac.uk>
;
;    * Modified to allow calls with only one argument
;

FUNCTION int_func, xin, fin, simple=simple
   IF N_PARAMS() EQ 1 THEN BEGIN
     f = xin
     x = FINDGEN(N_ELEMENTS(f))
   ENDIF ELSE BEGIN
     f = fin
     x = xin
   ENDELSE
   n = N_ELEMENTS(f)

   g = FLTARR(n)

   IF KEYWORD_SET(simple) THEN BEGIN
     ; Just use trapezium rule
     
     g[0] = 0.0
     FOR i=1, n-1 DO BEGIN
       g[i] = g[i-1] + 0.5*(x[i] - x[i-1])*(f[i] + f[i-1])
     ENDFOR
   ENDIF ELSE BEGIN
     
     n2 = FIX(n/2)
     
     g[0] = 0.0
     FOR i=n2, n-1 DO BEGIN
       g[i] = INT_TABULATED(x[0:i], f[0:i])
     ENDFOR
     
     FOR i=1, n2-1 DO BEGIN
       g[i] = g[n-1] - INT_TABULATED(x[i:*], f[i:*])
     ENDFOR
   ENDELSE
   
   RETURN, g
END
