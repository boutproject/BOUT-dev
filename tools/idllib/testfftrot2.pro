pro testfftrot,theta = theta
  

  
  set_plot,'ps'
  filename = 'testfft.ps'
  device,filename=filename,/color,YSIZE=25,YOFFSET =1,landscape =0
  print,"printing to : ",filename
  !p.multi = [0,1,2]
  loadct,40
  !x.thick=2
  !y.thick=2
  !p.thick=2
  !p.charthick=2
  !p.charsize =1

  TWOPI = 2.0 * !PI
  Lx = 200.0
  Ly = 32.0
  Nx = 32
  Ny = 100
  dx = Lx/Nx
  dy = Ly/Ny
 
  x=indgen(Nx)*dx
  y= indgen(Ny)*dy

  xx = rebin(indgen(Nx)*dx,Nx,Ny)
  yy = transpose(rebin(indgen(Ny)*dy,Ny,Nx))

  ii = rebin(indgen(Nx),Nx,Ny)
  jj = transpose(rebin(indgen(Ny),Ny,Nx))

  jj2 = jj+ii*sin(2)
  ii2 = ii+(jj)*1
  ii2 = (ii+jj*1/5.) mod Nx
  ;ii3 = (jj/sin(30)) mod Nx
  
  theta = Nx/Ny * !pi/4
  
  ii3 = (ii+ (jj* tan(theta))) mod Nx
  ii4 = (ii- (jj* tan(theta))) mod Nx
  jj3 = (jj+ (ii * tan(!pi/2 - theta))) mod Ny

  zeta = yy+xx*sin(1)

  zeta2 = zeta-xx*sin(1)

  image = sin(.01*x*!pi)#(sin(1.0/16.0*!pi*y))
;+.1*cos(.2*!pi*y ));+ (randomu(4,Nx,Ny)-.5)/2.0
  
  image2 = interpolate(image,ii3,jj,cubic = -.5)
  image3 = interpolate(image,ii,jj3,cubic=-.5)
  
  help,image2,image,ii2

  fft1 = fft(image)
  
  kx = 1/Lx * indgen(Nx)
  ky = 1/Ly * indgen(Ny)

  !p.multi = [0,2,3]
  contour,image,ii,jj,/fill,NLEVELS = 50,title = 'grid space'
  oplot,ii,jj,psym =3

  contour,image,ii,jj2,/fill,NLEVELS = 50,title = 'grid space'
  oplot,ii,jj2,psym =3
  
  contour,image,ii,jj,/fill,NLEVELS = 50,title = 'new indecies'
  oplot,ii,jj,psym =3
  oplot,ii2,jj,psym =3,color = 75

  contour,image2,title ='image 2',/fill,NLEVELS = 50
  
  contour,image3,title ='image 3',/fill,NLEVELS = 50
  
  contour,image,x,y,/fill,NLEVELS = 50,title = 'original'
  contour,image,x,zeta,/fill,NLEVELS = 50,title = 'skewed'
  contour,image,x,zeta2,/fill,NLEVELS = 50,title = 'unskewed'

  contour,rot(image,20),rot(xx,20),rot(zeta2,20),/fill,NLEVELS = 50,title = 'rotation?'
   contour,rot(image,0),rot(xx,45),rot(zeta2,45),/fill,NLEVELS = 50,title = 'rotation?'
  if not keyword_set(theta) then theta = 0
 
   Nxx = round(2*Nx);assume that Nx Ny are even
   Nyy = round(2*Ny)
  
   tile_image = dblarr(Nxx,Nyy)

   tile_image[Nx/2:Nxx-1-Nx/2,Ny/2:Nyy-1-Ny/2] = image2
   tile_image = tile_image+shift(tile_image,Nx,0)
   tile_image = tile_image+shift(tile_image,0,Ny)
   contour,tile_image,/fill,NLEVELS = 50


   ii = rebin(indgen(Nxx),Nxx,Nyy)
   jj = transpose(rebin(indgen(Nyy),Nyy,Nxx))

   jj2 = jj+ii*sin(2)
   ii2 = ii+(yy)*1
   zeta = yy+xx*sin(1)
 
  ;; sx = 10
  ;; sy = 10
  ;; dx= dx/sx
  ;; dy = dy/sy
  ;; dxx = sqrt((dx*cos(theta*!pi/180))^2+(dy*sin(theta*!pi/180))^2)
  ;; dyy = sqrt((dx*sin(theta*!pi/180))^2+(dy*cos(theta*!pi/180))^2)

  ;; MX = sx*2*Nx
  ;; MY = sy*2*Ny

  ;; xx = indgen(MX)*dxx ;twice the size of the earlier domain
  ;; yy = indgen(MY)*dyy ;""


  ;; ;print,"total(pad_Fft): ",total(pad_fft)
  ;; ;pad_image =  congrid(fft(pad_fft,/inverse),MX,MY,/interp)
  ;; pad_image = congrid(tile_image,MX,MY,/interp)
  ;; image3D = rebin(pad_image,MX,MY,2)
  ;; ;clip_image = pad_image[MX/4:MX-1-MX/4,MY/4:MY-1-MY/4]

  ;; ;debug
  

  ;; new_image = (rot(pad_image,theta,missing = 0))[MX/4:MX-1-MX/4,MY/4:MY-1-MY/4]
  ;; new_image3d = (transform_volume(image3D,rotation=[0,0,theta]))[MX/4:MX-1-MX/4,MY/4:MY-1-MY/4,*]
  ;; ;new_image = clip_image


  ;; xx = xx[0:MX/2-1]
  ;; yy = yy[0:MY/2-1]
  ;; Lxx = MX/2 * dxx
  ;; Lyy = MY/2 * dyy

  ;; contour,new_image,xx,yy,/fill,NLEVELS = 10,title="new image" ;;simple rotation
  ;; help,new_image3d[*,*,0],xx,yy
  ;; contour,reform(new_image3d[*,*,0]),xx,yy,/fill,NLEVELS = 10,title="new image 3d"
  
  ;; kxx =  1/(Lxx) * indgen(MX/2)
  ;; kyy = 1/(Lyy) * indgen(MY/2)
  
  ;; print,MX,MY,Lxx,Lx
  ;; print,"total(pad_image): ",total(pad_image),pad_image[1,25]
  ;; print,"total(new_image): ",total(new_image)
  ;; fft_new  = fft(new_image,/double)
  ;; ;pad_fft2 = (rot(pad_fft,3))[50:149,50:149]

  ;; contour,(alog(fft1*conj(fft1)))[1:Nx/2,1:Ny/2],kx[1:Nx/2],ky[1:Ny/2],/fill,nlevels = 50,/xlog,/ylog

  ;; contour,(alog(fft_new*conj(fft_new)))[1:MX/2-1,1:MY/2-1],kxx[1:MX/2-1],kyy[1:MY/2-1],/fill,nlevels = 50,/xlog,/ylog

  ;; contour,(alog(fft_new*conj(fft_new)))[1:Nx/2,1:Ny/2],kxx[1:Nx/2],kyy[1:Ny/2],/fill,nlevels = 50,/xlog,/ylog
  ;; contour,(alog(fft1*conj(fft1)))[1:Nx/4,1:Ny/4],kx[1:Nx/4],ky[1:Ny/4],/overplot,/xlog,/ylog,nlevels = 1

  ;; plot,kxx[1:Nx/2],kx[1:Nx/2]

  

  
  
  ;; ;; plot,kyy[0:30],(total(fft_new*conj(fft_new),1))[0:30],/YLOG,psym = 2
  ;; ;; oplot,ky[0:30],(total(fft1*conj(fft1),1))[0:30],psym = 4,color = 70
 
 
  
  ;; plot,ky[0:Ny/2],(total(fft1*conj(fft1),1))[0:Ny/2],/YLOG,psym = 4,color = 70,TITLE="Ky spectrum"
  ;; oplot,kyy,(total(fft_new*conj(fft_new),1)),psym = 2
  ;; oplot,kxx,(total(fft_new*conj(fft_new),2)),color = 10,linestyle = 2
  

  ;; plot,ky[0:Ny/2],((((total(fft1*conj(fft1),1))))/((total(fft_new*conj(fft_new),1))))[0:Ny/2],psym=2,/YLOG

   
  ;; plot,ky[0:Ny/2],(((total(fft1*conj(fft1),1))/Nx)[0:Ny/2])/((total(fft_new*conj(fft_new),1))/(Mx/2.)),/YLOG

  ;; plot,kx[0:Nx/2],(((total(fft1*conj(fft1),2))/Ny)[0:Nx/2])/((total(fft_new*conj(fft_new),2))/(My/2.)),/YLOG

  ;; plot,kx[0:30],(total(fft1*conj(fft1),2))[0:30],/YLOG,psym = 4,color = 70,TITLE="Kx spectrum"
  ;; oplot,kxx,(total(fft_new*conj(fft_new),2)),color = 10,linestyle = 2
  ;; oplot,kyy,(total(fft_new*conj(fft_new),1))


  ;; ;; plot,kyy[0:30],(total(fft_new*conj(fft_new),1))[0:30],/YLOG
  ;; ;; oplot,ky,(total(fft1*conj(fft1),1)),psym = 4,color = 70
  
 
  ;; ;; plot,kxx[0:30],(total(fft_new*conj(fft_new),2))[0:30],/YLOG
  ;; ;; oplot,kx,(total(fft1*conj(fft1),2)),psym = 4,color = 70
  
  ;; ;; plot,kxx,total(fft_new*conj(fft_new),2),/YLOG,psym = 2
  ;; ;; plot,kxx[0:30],(total(fft_new*conj(fft_new),2))[0:30],/YLOG,psym = 2
  ;; ;; oplot,kx[0:30],(total(fft1*conj(fft1),2))[0:30],psym = 4,color = 70

  ;; ;; contour,(rot(pad_image,15))[50:149,50:149]
  ;; ;; contour,fft(fft3,/inverse)
 
  device,/close

  pdf_save_name = strcompress(strjoin(['testfft','.pdf']),/remove_all)
  spawn,strjoin(["ps2pdf ",filename, " ",pdf_save_name])
end
