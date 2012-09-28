; growth_spec.pro
;
; ***************
; PRO growth_spec
; ***************
;
; Calculate growth rate and frequency for the full spectrum
; of each kz-Fourier mode for each point in x,y,k,t
;
; Input: 
;   f = data array: 1d-4d
;   t = optional time array
;   boolean average  = perform an average over x and y 
;   boolean xaverage = perform an average over x
;   boolean yaverage = perform an average over y
;
; Output:
;   spectrum = fourier spectrum in z-coordinate (assumed to be 3rd dimension)
;   growthk  = growth rate
;   freqk    = rotation frequency in z-direction
;   lambdak  = growthk + I*freqk 
;
;
; Created:  2011/08/09  I. Joseph
;
; Modified:
;
;   2012/03/28 I. Joseph: Added header
;                         Added multi-dimensional & averaging capabillity
;

pro growth_spec, fz, time=time, spectrum=fk, lambdak=lambdak, growthk=growthk, freqk=freqk, ac=ac, average=average, yaverage=yaverage, xaverage=xaverage

  d = size(fz)

  avg = keyword_set(average) or keyword_set(xaverage) or keyword_set(yaverage)  
  avg_x = keyword_set(average) or keyword_set(xaverage)
  avg_y = keyword_set(average) or keyword_set(yaverage)

  if d[0] gt 2 and avg then begin
    case d[0] of
      3: begin
           favg = mean(fz,dim=1)
           growth_spec, favg, time=t, spectrum=fk, lambdak=lambdak, growthk=growthk, freqk=freqk, ac=ac
         end
      4: begin
           if avgy then favg = mean(fz,dim=2)
           if avgx then favg = mean(favg,dim=1)
	   growth_spec, favg, time=t, spectrum=fk, lambdak=lambdak, growthk=growthk, freqk=freqk, ac=ac
         end
        default: print, 'ERROR: growth_spec requires an fz that is 1d-4d'
    endcase  
  endif else begin

    fk = fltarr(size(fz,/dimensions))
    lambdak = fltarr(size(fz,/dimensions))
    
    if not keyword_set(time) then time = findgen(d[d[0]])

    case d[0] of
      1: lambdak = deriv(time,alog(fz))
      2: begin
           fk = fft(fz,dim=1)
	   for k=0, d[1]-1 do $
             lambdak[k,*] = deriv(time,alog(reform(fk(k,*))) )
         end     
      3: begin
           fk = fft(fz,dim=2)
	   for j=0, d[1]-1 do $
	   for k=0, d[2]-1 do $
             lambdak[j,k,*] = deriv(time,alog(reform(fk(j,k,*))) )
         end 
      4: begin
           fk = fft(fz,dim=3)
	   for i=0, d[1]-1 do $
	   for j=0, d[2]-1 do $
	   for k=0, d[3]-1 do $
             lambdak[i,j,k,*] = deriv(time,alog(reform(fk(i,j,k,*))) )
         end 
      default: print, 'ERROR: growth_spec requires an fz that is 1d-4d'
    endcase

    growthk=real_part(lambdak)
    freqk=imaginary(lambdak)
  endelse
end



