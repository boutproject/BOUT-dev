function fftrot,input,tan_phi,rt=rt,spectrum=spectrum,view = view

;given some r and t this function will return a grid where one of the
;dimension corresponds to the by-normal direction
;/spectrum keyword only returns the spectral information
  

  if keyword_set(rt) then begin
     r= rt[0]
     t= rt[1]
     ;print,"r: ",r," t: ",t
  endif
  ;print,"pow"
 
  cot_phi = 1./tan_phi
  sin_phi = sqrt(1.0/(1.0+tan_phi^2))
  cos_phi = sqrt(1.0/(1.0+cot_phi^2))

  ;help,data,input

                                ;one can simply pass a slice in r,t,
                                ;or a full 4D simuation array, if the
                                ;later choice is made one can futher
  sizes = size(*(input.data))         ;chose to try to rotate it for all t
                                ; and r or specify r and t to save memory
  dim = sizes[0]
 
;size info should be drawn from the d str
  Nr = (**input.d).Nx
  Ny = (**input.d).Ny
  N_zeta = (**input.d).Nz -1 
  Nt = (**input.d).Nt




if dim NE 2 then begin ;if the user is providign something other than a 2d image
   Nt = sizes[4]
   if keyword_set(rt) then begin ;but is also providign r and t coordinates
     ;image = temporary(image[r,*,*,t])
     image = ((*(input.data)))[r,*,*,t] ;pull up the image like this
  endif else begin
     print,'loading all data'
     image = (*(input.data)) ;otherwise pull up all the data in the image
  endelse
endif else begin ;otherwise whatever the user provided is a 2d image 
   image = (*(input.data)) 
endelse

 ; loadct,40
  ;help,input,/str
  ;; help,input.L
Lz = input.L.zeta
Ly = input.L.b


;print,"Lz: ",Lz
  ;; print,N_zeta,Ny
 
  
TWOPI = 2.0 * !PI

dz = Lz/N_zeta
dy = Ly/Ny


z = indgen(N_zeta)*dz
y = indgen(Ny)*dy
z = z

  
  fft1 = fft(image)
  
  kz = TWOPI/Lz * indgen(N_zeta)
  ;print,kz
  ky = TWOPI/Ly * indgen(Ny)

  ;N_turns =  max((**input.d).Zxy)*
 
  N_turns = round(max(((**input.d).Zxy*(**input.d).Btxy)/(2*!PI*(**input.d).Bpxy)))

  
  N_zeta2 = round(N_turns*N_zeta)     ;assume that N_zeta Ny are even
  Nyy = round(Ny)
  
  ;help,Nyy,N_zeta2,N_turns
;we could make an array of pointers that point back to the original image
  ;tile_image = dblarr(N_zetaz,Nyy,Nr,Nt)
  if keyword_set(rt)  then  begin
   
     tile_image = dblarr(Nyy,N_zeta2)
     tile_image[*,0:N_zeta-1] = image
     ;tile_image = tile_image+shift(tile_image,Ny,0)
     for i_turn = 0,N_turns-1 do begin
        ;help,tile_image,image
        tile_image = tile_image+shift(reform(tile_image),0,N_zeta)
     endfor
     
     ;contour,tile_image,nlevels = 50 ,/fill
    

  endif else begin
     print,'rotating for all time steps at once - may cause memory issues'
     
     tile_image = dblarr(Nr,Nyy,N_zeta2,Nt)
     tile_image[*,Ny/2:Nyy-1-Ny/2,N_zeta/2:N_zeta2-1-N_zeta/2,*] = image
     tile_image = tile_image+shift(tile_image,0,Ny,0,0)
     tile_image = tile_image+shift(tile_image,0,0,N_zeta,0)
  endelse
  


  phi = (180.0/!pi) * atan(tan_phi)
  print,phi

  sz = 4
  sy = 4
  dz= dz/sz
  dy = dy/sy
     
  MZ = sz*N_turns*N_zeta
  MY = sy*Ny

  ii = rebin(indgen(MY),MY,MZ)
  jj = transpose(rebin(indgen(MZ),MZ,MY))


  


  if keyword_set(rt) then begin
     new_image = dblarr(MY,MZ)
  endif else new_image = dblarr(Nr,MY,MZ,Nt)


 
  if keyword_set(rt) then begin
     r0 = r 
     r1 = r
     Nr = 1
  endif else begin 
     r0 = 0
     r1 = Nr-1
  endelse
                                ;prepare the structure that will be output
  rotated = {image:ptr_new(),zz:ptrarr(Nr),yy:ptrarr(Nr),$
             kzz:ptrarr(Nr),kyy:ptrarr(Nr),Lzz:ptrarr(Nr),Lyy:ptrarr(Nr)} 

  r_i =0
  sin_theta = max((**input.d).BPXY,dimension=2)/max((**input.d).BXY,dimension=2)
  print,"sin_theta: ",sin_theta
  for r = r0,r1 do begin        ;loop over all relevant r's - may be just one index value
     ;cos_phi[r]=1
     ;sin_phi[r]=0

     ii2 = (ii- (jj* tan_phi[r])) mod MY
     print,'rotating . '
     ;dzz = sqrt((dz*cos_phi[r])^2+(dy*sin_phi[r])^2)
     ;dyy = sqrt((dz*sin_phi[r])^2+(dy*cos_phi[r])^2)
     dyy = dy
     ;dzz = sqrt((dz*cos_phi[r])^2+(dy*sin_phi[r])^2)
     dzz = dz*sin_theta[r]; wrong but whatever for now
     ;figure out grid related quantities for a given r
     
     
     zz = indgen(MZ)*dzz        ;twice the size of the earlier domain
     yy = indgen(MY)*dyy        ;""
     ;zz = zz[0:MZ/2-1]
    ; yy = yy[0:MY/2-1]
     Lzz = MZ * dzz
     Lyy = MY * dyy
     kzz =  TWOPI/(Lzz) * indgen(MZ)
     kyy = TWOPI/(Lyy) * indgen(MY)
     
     print,max(kzz),max(kyy),Lzz/N_turns,Lyy

     rotated.zz[r_i]=ptr_new(zz)
     (rotated.yy[r_i]) = ptr_new(yy)
     (rotated.kzz[r_i]) = ptr_new(kzz)
     (rotated.kyy[r_i]) = ptr_new(kyy)
     (rotated.Lzz[r_i]) = ptr_new(Lzz)
     (rotated.Lyy[r_i]) = ptr_new(Lyy)

  
     if keyword_set(rt) then begin
        ;given some fixed t
       
        
        pad_image = congrid(tile_image,MY,MZ,/interp)
        ;help,pad_image,ii2,jj
        new_image = interpolate(pad_image,ii,jj,cubic = -.5)
                                ;new_image=(rot(pad_image,-phi[r],missing
                                ;= 0))[MY/4:MY-1-MY/4,MZ/4:MZ-1-MZ/4]
        ;skew the up-sampled image
        
     endif else begin
        pad_image = congrid(reform(tile_image[r,*,*,*]),MY,MZ,Nt,/interp) 
        new_image[r,*,*,*]=(transform_volume(pad_image,rotation=[0,0,-phi[r]],missing = 0))[MY/4:MY-1-MY/4,MZ/4:MZ-1-MZ/4,*]
     endelse
    
     r_i++
  endfor
    

  rotated.image = ptr_new(new_image)

  print,'about to return this stuff'
   return,rotated
 
 
end
