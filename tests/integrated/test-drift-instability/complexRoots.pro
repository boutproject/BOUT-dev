function complexRoots, s, t, pr, pi, order=order
;
;
; s is a parameter, s>0
;----------------------

 t=( sqrt(s^4 + 16*s^2) - s^2 )/2.

 pr=sqrt(t)
 pi=2*s/pr


 if not keyword_set(ORDER) then order=1

 case order of
  1: begin 
      w=complex(pr,pi-s)/2
     end

  2: begin
      w=complex(-pr,-pi-s)/2
     end

 endcase

;
;
return,w
end

pro getComplexRoots, s, w1, w2, plot=plot, linear=linear

 if keyword_set(LINEAR) then begin
    s=1e-3 + 1e2*findgen(1001)/1000. 
 endif else begin
    s=1e-3*10^(6*findgen(100)/101)
 endelse

 w1=complexRoots(s, order=1)
 w2=complexRoots(s, order=2)

 if keyword_set(PLOT) then begin 
  plot, s, float(w1), yr=[-2,2], /xst, tit='Complex roots' 
  oplot, s, float(w2)
  oplot, s, imaginary(w1), lin=2 
  oplot, s, imaginary(w2), lin=2 &  oplot, s, s*0
 endif

end
