; a small procedure to remind me how to find the phase between two
; signal


pro wrap_phase

  ni = collect(var="Ni",path="data_1")

  f

end

pro find_phase_2D
  !p.multi = [0,1,1]
  phi = indgen(1000,/double)/100. ; a 10 second period, .01 second grid, or more generally a phase 

  x = indgen(32,/double); a 32 x32 spatial domain
  y = indgen(32,/double)


  ; the longest wavelength that can be represented by this 
  ; setup is L = 32, the shortest L = 4
  x1 = sin((1./16.)*(2.*!pi) * x) ; wavelength here is L_x = 16
  x2 = sin((1./16.)*(2.*!pi) * (x-2.)) ;another one shifted over by 2 space units,or 2*pi*2/16 =  pi/4 
  
  y1 = sin((1./8.)*(2.*!pi) * x) ; wavelength here is L_x = 8
  y2 = sin((1./8.)*(2.*!pi) * (x-2.)) ;here there is a shift of pi/2 

  ni1 = x1#y1
  ni2 = x2#y2

  ni1_fft = fft(ni1)
  ni2_fft = fft(ni2)

  
 

  pow1 = ni1_fft*conj(ni1_fft)
  
  pow2 = ni2_fft*conj(ni2_fft)
  
  cross_trans = conj(ni1_fft)*ni2_fft ;x1 leads x2
  
  cross_trans = cross_trans[0:16,0:16]
  phi = ATAN(cross_trans, /PHASE)
  amp = (abs(cross_trans))


  modes= array_indices(amp,where(max(amp)EQ amp))
 
  plot,[0,16],[0,16],psym =3
  oplot,[modes[0]],[modes[1]],psym = 2
  print,phi(where(max(amp)EQ amp))* (180./!pi); this is the net phase diff, scalar
  ;we can also find the vector of propagation
  dx = 1
  dy = 1
  dt  = 1
  v0 = [dx/dt,dy/dt]

  v_ph = v0*array_indices(amp,where(max(amp)EQ amp))
  print,v_ph

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


function find_phase_2D,slice1,slice2,dx = dx, dt= dt, dim1 = dim1, dim2 = dim2,show= show
  ;!p.multi = [0,1,1]
  
  loadct,40
  x = indgen(32,/double); a 32 x32 spatial domain
  y = indgen(32,/double)

  if not KEYWORD_SET(dx) then dx =1 
  if not KEYWORD_SET(dx) then dt =1
  if not KEYWORD_SET(dim1) then dim1 =1
  if not KEYWORD_SET(dim2) then dim2 =2

  N = (size(slice1))[dim1]
  M = (size(slice1))[dim2]
  print,"N x M:",N,"x",M


  ni1_fft = fft(slice1)
  ni2_fft = fft(slice2)

  
 

  pow1 = ni1_fft*conj(ni1_fft)
  
  pow2 = ni2_fft*conj(ni2_fft)
  
  cross_trans = conj(ni1_fft)*ni2_fft ;x1 leads x2
  
  cross_trans = cross_trans[0:N/2-1,0:M/2-1]
  phi = ATAN(cross_trans, /PHASE)
  amp = (abs(cross_trans))


  modes= array_indices(amp,where(max(amp)EQ amp))
 
  plot,[0,N/2-1],[0,M/2-1],psym =3
  ;contour,amp
  oplot,[modes[0]],[modes[1]],psym = 2
  print,phi(where(max(amp)EQ amp))* (180./!pi); this is the net phase diff, scalar
  ;we can also find the vector of propagation
  dx = 1
  dy = 1
  dt  = 1
  v0 = [dx/dt,dy/dt]

  v_ph = v0*array_indices(amp,where(max(amp)EQ amp))
  print,v_ph
  output = {modes:modes,v_ph:v_ph,amp:amp}
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



; if the dimension vector is not given simply return fft acroos all
; dims

; if the dimension vector IS given return fft AND 
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

  


  

  for i = 0, N_ELEMENTS(dimension)-1 do begin
     print,dimension[i]
     
     fft_input = (fft(fft_input,dimension = dimension[i]))

                                ;trim the excess, these does not seem
                                ;to be an easy way to do this
                                ;dynamically in idl
     
  endfor
  
  ;return, fft_input
                               ; sloppy non-general solution
  fft_input = fft_input[*,0:N[1]/2.,0:N[2]/2.,*]
  pow = fft_input * conj(fft_input)

  trash = max(total(abs(input),4),max_i)

  max_i = array_indices(total(input,4),max_i)

;by default the code search for the place with larges turbulence and
;perform wavelet analysis there
  if KEYWORD_SET(wavelet) then begin

     
     data =  reform(input[max_i[0],max_i[1],max_i[2],*])
     plot,data
     oplot,smooth(data,20),color = 75
     alt_data = data-smooth(data,20)

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

  pow = fft_input * conj(fft_input)
  
  ;for now assume that time is the last dimension 
  shift_param[Ndim-1] = 1.0

  cross_pow = fft_input * conj(shift(fft_input,shift_param))
  amp = abs(cross_pow)


  ;doing this is in a sloppy nongeneral way

  ;maximum across one dimension
  ;max_vals = max(amp,i,dimension=2)
  
  ;maximum across the 2nd dimension
 ; max_vals = max(max_vals,i,dimension=2)


                                ;nice way to get maximum values, but
                                ;does nto return indecies, make that
                                ;work later
  
  max_vals = cmapply('max',amp,[2,3]) ;see if we can we can provide a user define function for this 
  
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
           val= array_indices(reform(amp[i,*,*,j]),where(max_vals[i,j] EQ amp[i,*,*,j]))
           harmonic[i,j].loc = array_indices(reform(amp[i,*,*,j]),where(max_vals[i,j] EQ amp[i,*,*,j]))
           harmonic[i,j].phi = ATAN(cross_pow[i,val[0],val[1],j], /PHASE)
                                ;find the peaking wavenumber at any
                                ;given point in space and time

                                ;compute the phase diff at those peaked
                                ;peaked wavenumber between adjecent
                                ;slices in time

           harmonic[i,j].omega =  harmonic[i,j].phi/dt
           harmonic[i,j].pow = max_vals[i,j]
           
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

  Nx = (size(temp))[0]
  Nt = (size(temp))[1]
  
  print,keyword_set(wavelet)

  
  if keyword_set(wavelet) then output = {fft:fft_input,pow:pow,k:harmonic,wave_res:wave_res}
  
  
  if not keyword_set(wavelet) then output = harmonic



  ;oplot,(output.k.omega)[3,*],thick = 5
  
  ;temp = rebin(temp,N[0],N[3],
  ;modes= array_indices(amp,where(temp EQ amp[i,*,*,j]))

  if keyword_set(show) then begin
     if keyword_set(wavelet) then begin
     contour,reform(output.wave_res.pow[20,*])##reform(output.pow[3,*,5,20]),indgen(17),periods,/fill,nlevels= 20,$
             ystyle=1,/ylog,xstyle =1
     contour,reform(output.wave_res.pow[20,*])##reform(output.pow[3,*,5,20]),indgen(17),periods,/overplot
     
     contour,reform(output.wave_res.pow[50,*])##reform(output.pow[3,*,5,50]),indgen(17),periods,/fill,nlevels= 20,$
             ystyle=1,/ylog,xstyle =1
     contour,reform(output.wave_res.pow[50,*])##reform(output.pow[3,*,5,50]),indgen(17),periods,/overplot
     endif
  

  
 

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
