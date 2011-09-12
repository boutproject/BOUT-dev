; Filter for Discrete Cosine Transform 
; 
; Author:   Ilon Joseph, LLNL
;
; Began:    2011/06/02
; Modified: 2011/06/02


FUNCTION dctfilter, fin, tol=tol

        IF NOT KEYWORD_SET(tol) THEN tol=1e-3

        maxf = max(abs(fin))
        print, tol, maxf, maxf*tol

        fout = fin
        fout(WHERE( abs(fin) LE tol*maxf )) = 0.0
        
        return, fout
END
