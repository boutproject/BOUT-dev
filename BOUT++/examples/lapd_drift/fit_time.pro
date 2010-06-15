; Fit the solution from BOUT with a sine wave and calculate gamma/omega
; Work in progress

function RMS_XYZ, d
;Takes the rms over all space coordinates
;----------------------------------------------------------

;-find what tags are available in structure d
name=['APAR_XYZT','AJPAR_XYZT','JPAR_XYZT','NI_XYZT','PHI_XYZT','RHO_XYZT','TE_XYZT','TI_XYZT','VI_XYZT']
name_count=N_ELEMENTS(name)
iname=intarr(name_count)-1 ;-index of given name in structure d

dnames=TAG_NAMES(d)
nn=n_elements(dnames)
 
 for j=0,nn-1 do begin
  for i=0,name_count-1 do begin
    if (dnames[j] eq name[i]) then iname[i]=j
  endfor
 endfor


 dum=dblarr(d.trange)
 rms=dblarr(d.trange)


drms={APAR_T:dum,AJPAR_T:dum,JPAR_T:dum,NI_T:dum,PHI_T:dum,$
RHO_T:dum,TE_T:dum,TI_T:dum,VI_T:dum}

   for i=0,name_count-1 do begin

      IF iname[i] ge 0 THEN BEGIN
          for it=0, d.trange-1 do begin
            rms[it] = SQRT(TOTAL(double(d.(iname(i))[*,*,*,it])^2))
          endfor
          drms.(i)=rms
      ENDIF else begin
          drms.(i)=dum
      ENDELSE

      if not(keyword_set(QUIET)) then print, name(i), '...done'
   endfor

if keyword_set(DEBUG) then STOP

return, drms
end
;====================================================================================
; Return function for LMFIT:  
function fitFunc, X, A
; Fit the signal with
; f(x) = (a[0]*sin(a[1]+a[2]*x)

   vsin = a[0]*sin(a[1]+a[2]*x)
   return,[ [vsin],                    $ ; f(x)
            [vsin/a[0]],               $ ; d f(x) / d a[0]
            [a[0]*cos(a[1]+a[2]*x)],   $ ; d f(x) / d a[1]
            [x*a[0]*cos(a[1]+a[2]*x)]  $ ; d f(x) / d a[2]
          ]  
end  
;====================================================================================
function lmfit_signal, X, Y, A

;print, 'Initial fit parameters: ', A

measure_errors = Y

ffit = LMFIT(X, Y, A, /DOUBLE, CONVERGENCE=CONVERGENCE, MEASURE_ERRORS=measure_errors, $  
   FUNCTION_NAME = 'fitFunc', TOL=1.e-12)  

;print, 'Convergence=', CONVERGENCE
;print, 'Fit parameters: ', A

return, ffit

end
;====================================================================================
pro get_guess, X, Y, A
  
nz=zero_cross(Y, pos=pos)

if (nz > 0.) then begin
; Periodic data found, calculate the initial guess parameters

;print, 'nz=',nz
;print, 'Number of zeros = ', nz
;print, 'Positions of zeros:', pos
NT = (X[pos(nz-1)]-X[pos(0)])*2./(nz-1)
;print, 'Period:', NT

A = fltarr(3)
A[0] = max(Y)  ; amplitude
A[2] = 2.*!PI/NT ; omega

if (Y[0]<0.) then begin
  A[1] = - A[2]*X[pos(0)]  ; phase
endif else begin
  A[1] = !PI - A[2]*X[pos(0)]  ; phase
endelse

;print, A, X[pos(0)]
;
;Ndata = n_elements(Y)
;guess = fltarr(Ndata)
;for i=0,Ndata-1 do begin
;  afit=fitFunc(X[i],A)
;  guess[i] = afit[0]
;endfor

;plot, X, guess, title='Initial guess'

endif else begin

; No periodic data found (no zeros)
A = fltarr(3)
A[2] = -1.

endelse

end
;====================================================================================
pro SmoothData, tdata

  nt =  n_elements(tdata)
  ; Assume several oscillation periods are present in 0..nt-1 time interval
  ; Average out non-oscillating part if present
  tdata = tdata - smooth(smooth(tdata, nt/10,/edge_truncate), nt/10,/edge_truncate)

end
;====================================================================================

pro Fit_Time, d, du, MANUAL=MANUAL, NT0=NT0, NTMAX=NTMAX, omre=omre, omim=omim, $
                     SILENT=SILENT, PSAVE=PSAVE, LABEL=LABEL, SMOOTH=SMOOTH

 nstep=1
 nt=n_elements(d.PHI_xyzt[0,0,0,*])
 nx=n_elements(d.PHI_xyzt[*,0,0,0])
 ny=n_elements(d.PHI_xyzt[0,*,0,0])
 nz=n_elements(d.PHI_xyzt[0,0,*,0])

 if not keyword_set(LABEL) then LABEL='PLOT' ; label for the plot
 print, 'LABEL: ' + LABEL
 if not keyword_set(nt0) then nt0=nt/6 ;-skip initial part of the curve
 if not keyword_set(ntmax) then ntmax=nt ; include all points to the last one
 maxVal=dblarr(ntmax-nt0)

 rms=RMS_XYZ(d)
 for i=nt0,ntmax-1 do maxVal[i-nt0]=rms.PHI_t[i]

 if keyword_set(MANUAL) then begin
   plot, d.t_array[nt0:ntmax-1]/d.wci, alog(maxVal), psym=4, syms=3, xtit='time, s', ytit='ln<Ni>',/yst, chars=1.5

   print, "Mark 2 points on straight line to measure the exponential growth rate"
   print, "Click point 1" & mark, x=x1, y=y1 
   print, "Click point 2" & mark,x=x2,y=y2 
   oplot, [x1,x2], [y1,y2], col=2
   gamma=(y2-y1)/(x2-x1)
 endif else begin
   xx=d.t_array[nt0:ntmax-1]/d.wci & yy=ALOG(maxVal)
   res=Poly_Fit(xx,yy,1)
   oplot, xx, res[0] + res[1]*xx, col=2
   gamma=res[1]
   WAIT, 0.7
 endelse

 print, "gamma/OmCI=", gamma/d.wci

!p.multi=[0,1,2,0,0]
 WINDOW, RETAIN = 2
WINDOW, 0, XSIZE=750, YSIZE=750


 PHI_norm=d.PHI_XYZT[*,*,*,nt0:ntmax-1]
 PHI_fft=complexarr(nx,ny,nz,ntmax-nt0)
 growth=fltarr(ntmax-nt0)

 
 for i=nt0,ntmax-1 do growth[i-nt0]=exp(gamma/d.wci*d.T_ARRAY[i])
 for i=nt0,ntmax-1 do PHI_norm[*,*,*,i-nt0]=d.PHI_XYZT[*,*,*,i]/growth[i-nt0]

 jx=nx/3  ; choose some point away from the edges and not in the middle (avoid nodes of the wave)
 jy=ny/3
 jz=nz/3

 plot, d.t_array[nt0:ntmax-1]/d.wci, alog(maxVal), psym=4, syms=3, xtit='time, s', ytit='ln<Phi>',/yst, chars=1.5

 vd = PHI_norm[jx,jy,jz,*]
 SmoothData, vd

 plot, d.T_ARRAY[nt0:ntmax-1], vd, title='Phi(t)/exp(t)'

 omega_3d = fltarr(nx,ny,nz)

 for jx=0, nx-1 do begin

   if not keyword_set(SILENT) then print, 'jx = ', jx, '/', nx-1

   for jy=0, ny-1 do begin
   for jz=0, nz-1 do begin

     if keyword_set(SMOOTH) then begin ; de-trend
       vtmp = reform(PHI_norm[jx,jy,jz,*])
       SmoothData, vtmp
       PHI_norm[jx,jy,jz,*] = vtmp[*]
     endif
     get_guess, d.T_ARRAY[nt0:ntmax-1], PHI_norm[jx,jy,jz,*], A
     if (A[2] > 0) then ffit=lmfit_signal(d.T_ARRAY[nt0:ntmax-1],PHI_norm[jx,jy,jz,*], A)
     omega_3d[jx,jy,jz] = A[2]
   endfor
   endfor
 endfor

 jx=nx/3
 jy=ny/3
 jz=nz/3

; omega = A[2]
; print, 'omega/OmCI=',omega


 hmax=2.* TOTAL(omega_3d[where(omega_3d gt 0.)])/(n_elements(where(omega_3d gt 0.)))
 hbins = 200
 h=histogram(omega_3d,MIN=0.0, MAX=hmax, NBINS=hbins)
 homega = findgen(hbins)*hmax/(hbins-1)
 maxv = max(h,jmax)
 h_dom = homega(jmax)


 print, 'BOUT: omega/OmCI   = (', h_dom, '  +', gamma/d.wci, ')' 

 ; Return the values found
 omre = h_dom
 omim = gamma/d.wci


; plot, homega, h

;---------------------------------------------------------------
 jx=nx/3
 jy=ny/3
 jz=nz/3
jt=ntmax-1


!p.multi=[0,3,2,0,0]


window, 0, XSIZE=750, YSIZE=750, RETAIN = 2, /Pixmap,/Free
pixID=!D.Window

 WSet, pixID


 ; Plot the RMS of the signal (including the growth exp)
 plot, d.t_array[nt0:ntmax-1]/d.wci, alog(maxVal), psym=4, syms=3, xtit='time, s', ytit='ln<Phi>',/yst, chars=1.5, title=LABEL

 ; Plot Phi(t) (excluding the growth exp)
 plot, d.T_ARRAY[nt0:ntmax-1], PHI_norm[jx,jy,jz,*], title='Phi(t)/exp(t)'


 tit=string(format='(%"Omega 3d distribution, max @ om=%9.2e ")', h_dom)
 xtit='omega/OmCI'
 plot, homega, h, xtit=xtit, tit=tit

 tit='BOUT: Axial profile'
 xtit="Z, m"	 
 ytit="Ni"
 yr = [ min( [d.PHI_XYZT[jx,*,jz,jt]*0.9, d.PHI_XYZT[jx,*,jz,jt]*1.1]  ), $
        max( [d.PHI_XYZT[jx,*,jz,jt]*0.9, d.PHI_XYZT[jx,*,jz,jt]*1.1]  ) ]
 plot, du.ZXY[jx,*], d.PHI_XYZT[jx,*,jz,jt], $
       xtit=xtit, ytit=ytit, tit=tit, chars=1.5, yrange=yr

 tit='BOUT: Phi radial profile'
 xtit="r/rho_s"	 
 ytit="Phi"
 plot, du.RXY[*,jy]/d.rho_s, d.PHI_XYZT[*,jy,jz,jt], $
       xtit=xtit, ytit=ytit, tit=tit, chars=1.5

 tit='BOUT: Ni radial profile'
 xtit="r/rho_s"	 
 ytit="Ni"
 plot, du.RXY[*,jy]/d.rho_s, d.PHI_XYZT[*,jy,jz,jt], $
       xtit=xtit, ytit=ytit, tit=tit, chars=1.5

 !p.multi=0


 if keyword_set(PSAVE) then begin
;-save the picture
   image=TVRD(/true)
   file=LABEL+'.jpg'
   write_jpeg, file, image,/true
 endif

 WSet
 Device, Copy=[0, 0, !D.X_Size, !D.Y_Size, 0, 0, pixID]

;---------------------------------------------------------------


end
