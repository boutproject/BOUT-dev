function RMSvalue, vec1d
;
; -get rms of a 1D signal
;------------------------

  nel=n_elements(vec1d)
  valav=TOTAL(vec1d)/nel 
  valrms=SQRT(TOTAL((vec1d-valav)^2)/nel)
  acvec=vec1d-valav

;
return, {valrms:valrms, valav:valav, acvec:acvec}
end



pro moment_xyzt, sig_xyzt, rms=rms, dc=dc, ac=ac
;
; Calculate moments of a 4d signal of (x,y,z,t), i.e,
; -RMS, i.e., a function of (x,y,t)
; -DC (average in z), i.e., a function of (x,y,t)
; -AC (DC subtracted out), i.e., a function of (x,y,z,t)
;-------------------------------------------------------------------

  ON_ERROR, 2                     ; return to caller
 

  d = SIZE(sig_xyzt, /dimensions)
  IF N_ELEMENTS(d) NE 4 THEN BEGIN
      PRINT, "Error: Variable must be 4D (x,y,z,t)"
      RETURN
  ENDIF
  
  siz=SIZE(sig_xyzt)
  rms=fltarr(siz(1),siz(2),siz(4))
  dc=fltarr(siz(1),siz(2),siz(4))
  if keyword_set(AC) then ac=fltarr(siz(1),siz(2),siz(3),siz(4))
  
  data = sig_xyzt
  IF NOT is_pow2(siz[3]) THEN BEGIN
      PRINT, "WARNING: Expecting a power of 2 in Z direction"
      
      IF is_pow2(siz[3]-1) AND (siz[3] GT 1) THEN BEGIN
          PRINT, " -> Truncating last point to get power of 2"
          data = data[*,*,0:(siz[3]-2),*]
          siz[3] = siz[3] - 1
      ENDIF
  ENDIF
  
  for ix=0, siz(1)-1 do begin
      for iy=0, siz(2)-1 do begin
          for it=0, siz(4)-1 do begin
              val=RMSvalue(REFORM(sig_xyzt[ix,iy,*,it]))
              
              rms[ix,iy,it]=val.valrms
              dc[ix,iy,it]=val.valav
              if keyword_set(AC) then ac[ix,iy,*,it]=[val.acvec,val.acvec[0]]
              
          endfor
      endfor
  endfor
;
;
end
