pro find_peaks, sig, tt, tMax=tMax
;
; -Locate the peaks of sig and find the period
;------------------------------------------------;

  ;-isolate domains containing peaks
  ii=where(sig gt 0.5*max(sig))
  nii=n_elements(ii)
  ;plot, tt, sig


  ;-locate first maximum     
   iFirst1=0
   for i=0,nii-2 do begin
     if ((ii[i+1]-ii[i]) gt 1) then begin 
       ;-this is a gap in sequence
       i0=i
       break
     endif
   endfor
   iLast1=ii[i0]
   ;print, 'iFirst1:iLast1=', iFirst1, iLast1
   ;oplot, tt[iFirst1:iLast1], sig[iFirst1:iLast1], psym=2


  ;-locate second maximum     
   iFirst2=ii[i0+1]
   for i=i0+1,nii-2 do begin
     if ((ii[i+1]-ii[i]) gt 1) then begin 
       ;-this is a gap in sequence
       break
     endif
   endfor
   iLast2=ii[i]
   ;print, 'iFirst2:iLast2=', iFirst2, iLast2
   ;plots, tt[iFirst2], sig[iFirst2], psym=4
   ;oplot, tt[iFirst2:iLast2], sig[iFirst2:iLast2], psym=2


   ;-second maximum defines the period
   indices=iFirst2 +  indgen(iLast2-iFirst2)
   iMax=where(sig[indices] eq max(sig[indices]))

   tMax=tt[indices[iMax]]
   sigMax=sig[indices[iMax]]
   print, 'tMax=', tMax

   plots, tMax, sigMax, psym=4, syms=4
;
;
;
end
