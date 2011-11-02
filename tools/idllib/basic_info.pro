function basic_info, data,rescale = rescale,name = name,time_x=time_x,x=x,wavelet=wavelet,rotate = rotate,L=L,info = d

  save_mem =1
;; from some experimenting it appear that I need to limit the range of
;; k's I consider to roughly a half of what I should ahve to

;help,*d,/str
help,(*d).BPXY
;print,"theta: ",theta
;; assumes 4D data

nx = N_ELEMENTS(data[*,0,0,0])
ny = N_ELEMENTS(data[0,*,0,0])
nz = N_ELEMENTS(data[0,0,*,0])
nt = N_ELEMENTS(data[0,0,0,*])
TWOPI = 2.0 * !PI

;; tan_theta =  mean(((*d).BPXY)[*,ny/2])/mean(((*d).BTXY)[*,ny/2])
;; tan_theta =  mean(((*d).BPXY)[*,ny/2])/mean(((*d).BTXY)[*,ny/2])

tan_theta =max((*d).BPXY,dimension=2)/max((*d).BTXY,dimension=2)
;tan_phi =
;(*d).DZ*TWOPI*(max((*d).Rxy,dimension=2)/max((*d).Zxy,dimension=2))*tan_theta*(*d).Ny/(*d).Nz
cos_theta = max((*d).BTXY,dimension=2)/max((*d).BXY,dimension=2)

tan_phi = (*d).DZ*TWOPI*max((*d).Rxy/(*d).ycoord,dimension=2)* cos_theta*(*d).Ny/(*d).Nz
cot_phi = 1./tan_phi
sin_phi = sqrt(1.0/(1.0+tan_phi^2))
cos_phi = sqrt(1.0/(1.0+cot_phi^2))


print,"d.Zeta: ",(*d).DZ
help,tan_theta,tan_phi
;; plot,tan_theta
;plot,atan(tan_phi)
;plot,atan(tan_theta)

;;Global DC
DC_g = cmapply('+',data,[1,2,3])/(nx*ny*nz) ;;global DC - inherently bad?

;;global amplitude 
amp_gt = cmapply('max',abs(data),[1,2,3]) 
;;not a bad way to do things so long as this is being use on
;;fluctuating quantities, moreover the DC comp really needs to be added
;;back onto the steady state piece during the simulation

if total(DC_g) eq 0 then tempDC =0 else tempDC = transpose(rebin(DC_g,nt,nx,ny,nz),[1,2,3,0])

;;DC_g offset amplitude
amp_go = cmapply('max',abs(data-tempDC),[1,2,3]) ;;subtract off the mean, this is the max for any given timestep

if total(amp_go) EQ 0 then begin
   amp_stat = 0.0
   rescale = 0.  ; turn off rescaling if there is no ac comp to norm w/ respect to
endif

four_D_amp = transpose(rebin(amp_go,nt,nx,ny,nz),[1,2,3,0])  ;;resize to agree with dimensions of data

;;all spectral analysis is handled by find_phase_ND
;;
; let's determine the angle we need to rotate the grid by to
; look at bynormal modes




if KEYWORD_SET(rescale) then begin

   ;;normalized global DC 
   DC_gn = DC_g/amp_go
   
  
   ;;normalized input data
   data_n = data/four_D_amp ;;let you see if you have a traveling wave withoug resetting the sc
   
   ave = {data_n:data_n,DC_gn:DC_gn,amp_gt:amp_gt,DC_g:DC_g,amp_go:amp_go}
   to_fft = ptr_new(data_n)
endif else begin
   outdata= data
   ave = {data:data,amp_gt:amp_gt,DC_g:DC_g,amp_go:amp_go}
   to_fft = ptr_new(data)
endelse
;if KEYWORD_SET(rescale) then data = data/transpose(rebin(amp_go,nt,ny,nz,nx))

if keyword_set(rotate) then begin
   print,'have to rotate this shit'

   input = {data:ptr_new(data),L:L,d:ptr_new(d)}
   
   N_turns = max(round(max(((*d).Zxy*(*d).Btxy)/(2*!PI*(*d).Bpxy))))
   print,"N_turns: ",N_turns
 
   if keyword_set(save_mem) then begin
                                ; its even clear that this method
                                ; saves memory, buts its not slower
                                ; . . so wtf
      ;rot_data = fftrot(input,tan_phi,rt=[0,0])
      ;rot_data =  replicate(rot_data,nx,nt)
      rot_data = {image:ptr_new(dblarr(nx,4*ny,4*N_turns*nz,nt)),zz:ptrarr(nx),yy:ptrarr(nx),$
                  kzz:ptrarr(nx),kyy:ptrarr(nx),Lzz:ptrarr(nx),Lyy:ptrarr(nx)} 
      for i = 0,nx-1 do begin
         print,i
         for j= 0,nt-1 do begin
            
            print,"j:",j
            ;rot_data[i,j] = fftrot(input,tan_phi,rt=[i,j]) 
            temp = fftrot(input,tan_phi,rt=[i,j])
            ;print,'back in basic_info.pro'
            ;; help,temp,/str
            ;; help,temp.yy
            ;; help,*(rot_data.image)
         
         
         
            rot_data.yy[i] = (temp.yy[0])
            rot_data.zz[i] = (temp.zz[0])
            rot_data.kyy[i] = (temp.kyy[0])
            rot_data.kzz[i] = (temp.kzz[0])
            rot_data.Lyy[i] = (temp.Lyy[0])
            rot_data.Lzz[i] = (temp.Lzz[0])
            ;help,rot_data.image,temp.image
            ;help,*(rot_data.image),*(temp.image)
            (*(rot_data.image))[i,*,*,j]=reform(temporary(*(temp.image)))
            
         endfor
      endfor
   endif else begin
      rot_data = fftrot(input,tan_phi) 
   endelse
                                
                                ;we really need to maniulate the
                                ;rotated data such that the out looks
                                ;the same regardless of how we
                                ;produced - this functions need to
                                ;           operate like a black box
                                ;           for procidures and
                                ;           routines that call it


   ;; contour,data[0,*,*,1],/fill,NLEVELS = 50
   ;; contour,(*(rot_data.image))[0,*,*,0],*(rot_data.yy[0]),*(rot_data.zz[0]),title='rotated image',/fill,NLEVELS = 100
   ;; contour,(*(rot_data.image))[3,*,*,50],*(rot_data.yy[3]),*(rot_data.zz[3]),title='rotated image',/fill,NLEVELS = 100
   
   ;; if keyword_set(save_mem) then begin
   ;;    print,'saving memory'
   ;;    ;; contour,*((rot_data[0,1]).image),*((rot_data[0,1]).yy[0]),*((rot_data[0,1]).zz[0]),title='rotated image',/fill,NLEVELS = 50
   ;;    ;; contour,*((rot_data[nx/2,1]).image),*((rot_data[nx/2,1]).yy[0]),*((rot_data[nx/2,1]).zz[0]),title='rotated image',/fill,NLEVELS = 100
   ;; endif else begin
   ;;    contour,(*(rot_data.image))[0,*,*,0],*(rot_data.yy[0]),*(rot_data.zz[0]),title='rotated image',/fill,NLEVELS = 100
   ;;    ;; contour,*((rot_data[nx/2,0]).image),*((rot_data[nx/2,0]).yy[0]),*((rot_data[nx/2,0]).zz[0]),title='rotated image',/fill,NLEVELS = 50
   ;; endelse
   
   ;plot,pad_data[round(nx/2),0,*,0]
   ;; device,/close
   
endif


fft_info = find_phase_ND(*to_fft,dimension =[2,3],wavelet=wavelet,show=show,x=x,t=time_x,rescale = rescale)



if not keyword_set(name) then name = ''


if keyword_set(rotate) then begin
   rot_fft_info = find_phase_ND(*(rot_data.image),dimension =[2,3],wavelet=wavelet,show=show,x=x,t=time_x,rescale = rescale)
help,fft_info,/str

   output = {ave:ave,fft:temporary(fft_info),nx:nx,ny:ny,nz:nz,nt:nt,rot:temporary(rot_data),rot_fft:temporary(rot_fft_info)}
endif else begin
   print,"outputing. . "
   output = {ave:ave,fft:temporary(fft_info),nx:nx,ny:ny,nz:nz,nt:nt}
endelse



;; pow_data= fft_info.pow
;; ky = (total(total(pow_data,3)/nz,1)/nx) ; averaged over kz and radial direction
;; kz = (total(total(pow_data,2)/ny,1)/nx) ; averaged over ky and radial

;; kz_t = total(pow_data,2)/ny ; averaged over ky only
;; r_pow = total(total(pow_data,2)/ny,2)/nz

;; pow_data= total(pow_data,1)/nx ; radial average


;; pow_data_full= abs(fft_data_full)^2
;; kz_full = (total(total(pow_data_full,2)/ny,1)/nx); averaged over ky
;; ky_full = (total(total(pow_data_full,3)/nz,1)/nx); averaged over ky

;; pow_data_full= total(abs(fft_data_full)^2,1)/nx ; radial average




 ;; if KEYWORD_SET(name) then output =  {nx:nx,ny:ny,nz:nz,nt:nt,DC:DC,amp:amp,amp_t:amp_t,$
;;                                       pow:{k:pow_data,ky:ky,kz:kz},fft:fft_data,data:outdata,name:name,status:1} $
;;  else  output =  {nx:nx,ny:ny,nz:nz,nt:nt,DC:DC,amp:amp,amp_t:amp_t,$
;;                   pow:{k:pow_data,ky:ky,kz:kz,kz_full:kz_full,ky_full:$
;;                        ky_full,full_pow:pow_data_full,r:r_pow,kz_t:kz_t},$
;;                   fft:fft_data,data:outdata,status:1}


;; output =  ptr_new({nx:nx,ny:ny,nz:nz,nt:nt,DC:DC,amp:amp,amp_t:amp_t,$
;;                                        pow:{k:pow_data,ky:ky,kz:kz},$
;;            harmonics:fft_info,data:outdata,name:name,status:1})

 return,output
end
