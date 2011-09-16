function find_phase_ND,input, dimension = dimension, x = x, t =t, short = short,$
                       nscales = nscales,wavelet=wavelet,show = show,rescale = rescale

 
  

  if keyword_set(show) then begin
     set_plot,'ps'
     save_name='phase_view.ps'
     device,filename=save_name,/color,YSIZE=25,YOFFSET =1,landscape = 0
     loadct,40
     !x.thick=3
     !y.thick=3
     !p.thick=3
     !p.charthick=3
     !p.charsize =2
     !p.multi = [0,1,2]
  endif
  
  if not KEYWORD_SET(dimension) then return,fft(input)

  ; from here on out we can assume that the dimension vector is set
  Ndim = (size(input))[0]
  N =  (size(input))[1:Ndim]

  if not keyword_set(x) then x = indgen(N[0],N[1],N[2],/double)
  if not keyword_set(t) then t = indgen(N[3],/double)

  dt = t[1] - t[0]

  fft_input = input
  shift_param =dblarr(Ndim)
  

                                ; the best way to do this is to take a
                                ; fft in y and z and apply wavelet
                                ; techiques in the temporal dimension
                                ; to see the corresponding omegas

  


  

  for i = 0, N_ELEMENTS(dimension)-1 do begin ; typically this is the 2nd and 3rd dimension
     print,dimension[i]
     
     fft_input = (fft(fft_input,dimension = dimension[i]))

                                ;trim the excess, these does not seem
                                ;to be an easy way to do this
                                ;dynamically in idl
     
  endfor
  
  ;return, fft_input
                               ; sloppy non-general solution

  f_ny = N[1]/2. -1
  f_nz = N[2]/2. -1
  fft_input = fft_input[*,0:f_ny,0:f_nz,*]
  pow = fft_input * conj(fft_input)

  trash = max(total(abs(input),4),max_i) ; maximum at any given point in time

  max_i = array_indices(total(input,4),max_i)

;by default the code searches for the place with largest amplitude and
;perform wavelet analysis there
  if KEYWORD_SET(wavelet) then begin

     
     data =  reform(input[max_i[0],max_i[1],max_i[2],*])
     plot,data
     smoothed_data  = smooth(data,N_ELEMENTS(data)/9.)
     oplot,smoothed_data,color = 75
     alt_data = data-smoothed_data

     data_norm = data/ smooth(cmapply('max',abs(input),[1,2,3]),3) ;; there has GOT to be a better way to create a time series for wavelet analysis
     pure_grow = reform(data)/max(input[*,*,*,0])-data_norm


     data = data_norm
     
     plot, cmapply('max',abs(input),[1,2,3])
     
     plot,data,TITLE = "Normalized response"
     
     plot, abs(reform(data)/max(input[*,*,*,0])),/ylog
     oplot,abs(pure_grow),color = 75
     oplot,abs(data),color = 135
     plot,(data),color = 135

                                ;k0 is the wavelet wavenumner, 6 by
                                ;default for morlet
     ;data= alt_data
     plot,alt_data

     k_0 = 6

    
     f_high = 1.0/(2*dt)
     f_low = 1.0/(N[3]*dt)
     fourier_factor = (4*!PI)/(k_0 + SQRT(2+k_0^2)) 
                                ;fourier_factor* scale = wavelength/waveperiod in the fourier sense
     
                                ; the f_high is arbitrary units, as long a correct dt is supplied below
     s0 = round((1.0/f_high)* (1/fourier_factor))
     
     print,s0
     
     s_max =round( (1.0/f_low)* (1/fourier_factor))
     
     if not KEYWORD_SET(n_scales) then n_scales = 20
     dj = 1.0/(n_scales*alog(2))* ( alog(s_max)-alog(s0))

;;   print,dj
     
     waveout=temporary(wavelet(data,dt,PERIOD=periods,FFT_THEOR=fft_theor,$
                               /pad,mother='morlet',param=k_0,s0=s0,dj=dj,j=n_scales,COI = COI,SIGNIF = SIGNIF))


     help,waveout

     help,signif,fft_theor
     
     plot,periods,signif
     plot,periods,fft_theor


     wave_res = {pow:waveout*conj(waveout),period:periods,omegas:(2*!PI)/periods}

     if keyword_set(show) then begin
        contour,waveout*conj(waveout),t,periods,/fill,nlevels= 20,/ylog,ystyle=1
        contour,waveout*conj(waveout),t,periods,/overplot
        
     
        PLOTS,t,coi,NOCLIP=0,linestyle = 2
        
        contour,waveout,t,periods,/fill,nlevels= 20,/ylog,ystyle =1 
        contour,waveout,t,periods,/overplot
        contour,waveout*conj(waveout),t,periods,/overplot,nlevels = 2,$
                C_THICK=[7,7],C_COLOR=[0,0] ;,C_LINESTYLE = 1
        
        PLOTS,t,coi,NOCLIP=0,linestyle = 2


        wave_res = {pow:waveout*conj(waveout),period:periods,omegas:(2*!PI)/periods}
        contour,waveout*conj(waveout),t,(2*!PI)/periods,/fill,nlevels= 20,/ylog,ystyle=1
        contour,waveout*conj(waveout),t,(2*!PI)/periods,/overplot
        
     
        PLOTS,t,(2*!PI)/coi,NOCLIP=0,linestyle = 2
        
        contour,waveout,t,(2*!PI)/periods,/fill,nlevels= 20,/ylog,ystyle =1 
        contour,waveout,t,(2*!PI)/periods,/overplot
        contour,waveout*conj(waveout),t,(2*!PI)/periods,/overplot,nlevels = 2,$
                C_THICK=[7,7],C_COLOR=[0,0] ;,C_LINESTYLE = 1
        
        PLOTS,t,(2*!PI)/coi,NOCLIP=0,linestyle = 2
     endif
  endif

  
  ;for now assume that time is the last dimension 
  shift_param[Ndim-1] = 1.0

  cross_pow = fft_input * conj(shift(fft_input,shift_param))
  pow = fft_input * conj(fft_input)

  pow_kz = cmapply('+',pow,2)/f_ny
  pow_ky = cmapply('+',pow,3)/f_nz

  pow_rt = cmapply('+',pow,[2,3])/(f_ny*f_nz)

  pow_kzt = cmapply('+',pow,[1,2])/(N[0]*f_ny)
  pow_kyt = cmapply('+',pow,[1,3])/(N[0]*f_nz)
  
  pow = {pow:pow,kz:pow_kz,ky:pow_ky,rt:pow_rt,cross:cross_pow,kyt:pow_kyt,kzt:pow_kzt}
  
  amp = abs(cross_pow)

  ;doing this is in a sloppy nongeneral way

  ;maximum across one dimension
  ;max_vals = max(amp,i,dimension=2)
  
  ;maximum across the 2nd dimension
 ; max_vals = max(max_vals,i,dimension=2)


                                ;nice way to get maximum values, but
                                ;does nto return indecies, make that
                                ;work later
  
  max_vals = cmapply('max',amp,[2,3]) ;see if we can we can provide a user defined function for this later
  
  ;we can find higher order modes .. but for now
  N_modes = 1
  

  ;max_vals = rebin(max_vals,N[0],N[3],N[1]/2+1,N[2]/2+1)
  ;max_vals = transpose(max_vals,[0,2,3,1])

  ;modes = {y_mode:lonarr(N[0],N[3]),z_mode:lonarr(N[0],N[3]),phi:dblarr(N[0],N[3])}

  harmonic = replicate({loc:lonarr(2),phi:0.0,omega:0.0,pow:0.0},N[0],N[3])
  

  ; write an array based version of this later, for now N_modes = 1
  for k = 0,N_modes-1 do begin
     for i = 0,N[0]-1 do begin
        for j = 0,N[3]-1 do begin

           if max_vals[i,j] NE 0 then begin
              val= array_indices(reform(amp[i,*,*,j]),where(max_vals[i,j] EQ amp[i,*,*,j]))
              harmonic[i,j].loc = val[*,0]
              harmonic[i,j].phi = ATAN(cross_pow[i,val[0],val[1],j], /PHASE)
                                ;find the peaking wavenumber at any
                                ;given point in space and time
              
                                ;compute the phase diff at those peaked
                                ;peaked wavenumber between adjecent
                                ;slices in time
              
              harmonic[i,j].omega =  harmonic[i,j].phi/dt
              harmonic[i,j].pow = max_vals[i,j]
              
           endif else begin
              val= [0,0]
              harmonic[i,j].loc = [0,0]
              harmonic[i,j].phi = 0
                                ;find the peaking wavenumber at any
                                ;given point in space and time
              
                                ;compute the phase diff at those peaked
                                ;peaked wavenumber between adjecent
                                ;slices in time
              
              harmonic[i,j].omega = 0
              harmonic[i,j].pow = 0
              
           endelse
                                ;assuming there is single omega
                                ;responsibel for this compute it
        endfor
     endfor
  endfor
  
                    ;now we have to loop over the non-fft
                                ;dimensions and identify the position
                                ;in space and time  . . ? have to thnk
                                ;mor abot this
  ;now we have a 

 ;;  Nx = (size(temp))[0]
;;   Nt = (size(temp))[1]
  
;;   print,keyword_set(wavelet)

  
  
 
  
  
  if not keyword_set(wavelet) then wave_res = 0.0

  
  
  output = {pow:pow,harmonics:harmonic,wave_res:wave_res,ny:f_ny,nz:f_nz}



  ;oplot,(output.k.omega)[3,*],thick = 5
  
  ;temp = rebin(temp,N[0],N[3],
  ;modes= array_indices(amp,where(temp EQ amp[i,*,*,j]))

  if keyword_set(show) then begin
    ;;  if keyword_set(wavelet) then begin
;;      contour,reform(output.wave_res.pow[20,*])##reform(output.pow[3,*,5,20]),indgen(17),periods,/fill,nlevels= 20,$
;;              ystyle=1,/ylog,xstyle =1
;;      contour,reform(output.wave_res.pow[20,*])##reform(output.pow[3,*,5,20]),indgen(17),periods,/overplot
     
;;      contour,reform(output.wave_res.pow[50,*])##reform(output.pow[3,*,5,50]),indgen(17),periods,/fill,nlevels= 20,$
;;              ystyle=1,/ylog,xstyle =1
;;      contour,reform(output.wave_res.pow[50,*])##reform(output.pow[3,*,5,50]),indgen(17),periods,/overplot
;;      endif
  

  
 

     device,/close
     set_plot,'X'
     !x.thick=1
     !y.thick=1
     !p.thick=1
     !p.charthick=1
     !p.charsize =1
     !p.multi = [0,1,1]
     pdf_save_name = 'phase_view.pdf'
     spawn,strjoin(["ps2pdf ",save_name, " ",pdf_save_name])
  endif
  
;;   cross_trans = conj(ni1_fft)*ni2_fft ;x1 leads x2
  
;;   cross_trans = cross_trans[0:N/2-1,0:M/2-1]
;;   phi = ATAN(cross_trans, /PHASE)
;;   amp = (abs(cross_trans))


;;   modes= array_indices(amp,where(max(amp)EQ amp))
 
;;   plot,[0,N/2-1],[0,M/2-1],psym =3
;;   ;contour,amp
;;   oplot,[modes[0]],[modes[1]],psym = 2
;;   print,phi(where(max(amp)EQ amp))* (180./!pi); this is the net phase diff, scalar
;;   ;we can also find the vector of propagation
;;   dx = 1
;;   dy = 1
;;   dt  = 1
;;   v0 = [dx/dt,dy/dt]

;;   v_ph = v0*array_indices(amp,where(max(amp)EQ amp))
;;   print,v_ph
;;   output = {modes:modes,v_ph:v_ph,amp:amp}
  return,output
  ;use one of these signals are a reference
 ;;  ref_sig_info = primary_harmonic(x1,t,5)

;;   plot,ref_sig_info.fscale,ref_sig_info.f_autocorr
  
;;   period = ref_sig_info.f0
;;   peaks = ref_sig_info.peaks
;;   F = ref_sig_info.Fscale
;;   mag = ref_sig_info.f_autocorr
;;   harmonic_index = ref_sig_info.harmonic_i
;;   print,harmonic_index
;;   print,'peaks: ',peaks
;;   print,'period: ',period
;;   print,'phase diff: ',phi[harmonic_index]*180./!pi

end
