; a small procedure to remind me how to find the phase between two
; signal

pro find_phase
  !p.multi = [0,1,2]
  t = indgen(1000,/double)/100. ; a 10 second period, .01 second grid
  x1 = sin(2.*(2.*!pi) * t) ; with a period of .5 second, or freq of 2Hz
  x2 = sin(2.*(2.*!pi) * (t-.125) ) ;another one shifted over by .125 seconds,or 90 deg,pi/2  

  plot,t[0:100],(x1)[0:100]
  oplot,t[0:100],x2[0:100],color = 75,linestyle =2 

  fft_x1 = fft(x1)
  fft_x2 = fft(x2)

  pow1 = fft_x1*conj(fft_x1)
  
  pow2 = fft_x2*conj(fft_x2)
  
  cross_trans = conj(fft_x1)*fft_x2 ;x1 leads x2
  N = N_ELEMENTS(cross_trans)

  phi = ATAN(cross_trans[0:N/2.], /PHASE)
  amp = (abs(cross_trans))

  ;use one of these signals are a reference
  ref_sig_info = primary_harmonic(x1,t,5)

  plot,ref_sig_info.fscale,ref_sig_info.f_autocorr
  
  period = ref_sig_info.f0
  peaks = ref_sig_info.peaks
  F = ref_sig_info.Fscale
  mag = ref_sig_info.f_autocorr
  harmonic_index = ref_sig_info.harmonic_i
  print,harmonic_index
  print,'peaks: ',peaks
  print,'period: ',period
  print,'phase diff: ',phi[harmonic_index]*180./!pi

end
