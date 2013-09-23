;Import the experimentally measured temperature profiles, inarr, into the
;grid file, g. The input profile should be interpolated to the grid
;before this code.(ver 0.1). This function is only for single null geometry

function smoothtwo, inarr, xmin=xmin, xmax=xmax, width=width
  ON_ERROR, 2
  
  nx=n_elements(inarr)
  p=inarr
  ptemp=p[xmin:xmax]
  pt5=smooth(ptemp,width)
  pt5=smooth(pt5,5)
  for x=xmin,xmax do begin
     p[x]=pt5[x-xmin]
  endfor
  return,p
end

pro Ni2Gridalls, fileni, filete, fileti, filename, smooth=smooth 
  
  ON_ERROR, 2

  if ( size(fileni,/type) ne 7 ) then begin 
     print, "Wrong ion density input file."
     return
  endif

  if ( size(filete,/type) ne 7 ) then begin 
     print, "Wrong electron temperature input file."
     return
  endif

  if ( size(fileti,/type) ne 7 ) then begin 
     print, "Wrong ion temperature input file."
     return
  endif

  if ( size(filename,/type) ne 7 ) then begin 
     print, "Wrong grid file."
     return
  endif  

  g=file_import(filename)
  rm_temp=max(g.rxy[g.nx-1,*],ypeak)
  psn=(g.psixy[*,ypeak]-g.psi_axis)/(g.psi_bndry-g.psi_axis)

  restore,fileni
;  ni=profp
  restore,filete
;  te=profp
  restore,fileti

  nx = g.nx
  ny = g.ny

  kb = 1.38e-23
  profte = fltarr(nx,ny)
  profni = fltarr(nx,ny)
  profti = fltarr(nx,ny)
  profp = fltarr(nx,ny)
  xsep = floor(g.ixseps1)

  profp = g.pressure
  tmp = where(te gt 0.01)
  nte = n_elements(tmp) 
  tet = te
  for i=nte, nx-1 do begin
     tet[i] = te[nte-1]
  endfor
;  tet=smoothtwo(tet,xmin=nte-20,xmax=nte+20,width=20)

  for i=0, xsep do begin
     for j=0, ny-1 do begin
        if ((j gt g.jyseps1_1) and (j le g.jyseps2_1)) or ((j gt g.jyseps1_2) and (j le g.jyseps2_2)) then begin
           profte[i,j] = tet[i]
           profni[i,j] = ni[i]
        endif else begin
           profte[i,j] = tet[xsep]
           profni[i,j] = ni[xsep]
        endelse
     endfor
  endfor

  for i=xsep, nx-1 do begin
     for j=0, ny-1 do begin
        profte[i,j] = tet[xsep]
        profni[i,j] = ni[xsep]
     endfor
  endfor

  profti = g.pressure/(16.02*profni)-profte

  if keyword_set(smooth) then begin
     xmin1=xsep-20
     xmax1=xsep+20
     for j=0,ny-1 do begin
        profp[*,j]=smoothtwo(profp[*,j],xmin=xmin1,xmax=xmax1,width=20)
        profni[*,j]=smoothtwo(profni[*,j],xmin=xmin1,xmax=xmax1,width=20)
        profte[*,j]=smoothtwo(profte[*,j],xmin=xmin1,xmax=xmax1,width=20)
        profti[*,j]=smoothtwo(profti[*,j],xmin=xmin1,xmax=xmax1,width=20)
     endfor
  endif

  for i=xsep, nx-1 do begin
     for j=0, g.jyseps1_1 do begin
        profte[i,j] = profte[i,g.jyseps1_1+1]
        profti[i,j] = profti[i,g.jyseps1_1+1]
        profni[i,j] = profni[i,g.jyseps1_1+1]
        profp[i,j] = profp[i,g.jyseps1_1+1]
     endfor
     for j = g.jyseps2_2+1, ny-1 do begin
        profte[i,j] = profte[i,g.jyseps2_2]
        profti[i,j] = profti[i,g.jyseps2_2]
        profni[i,j] = profni[i,g.jyseps2_2]
        profp[i,j] = profp[i,g.jyseps2_2]
     endfor
  endfor

  for i=0, xsep-1 do begin
     for j=0, g.jyseps1_1 do begin
        profte[i,j] = profte[xsep,j]
        profti[i,j] = profti[xsep,j]
        profni[i,j] = profni[xsep,j]
        profp[i,j] = profp[xsep,j]
     endfor
     for j = g.jyseps2_2+1, ny-1 do begin
        profte[i,j] = profte[xsep,j]
        profti[i,j] = profti[xsep,j]
        profni[i,j] = profni[xsep,j]
        profp[i,j] = profp[xsep,j]
     endfor
  endfor  

;  f=file_export(filename, gnew)
  handle = file_open(filename, /write)
  s = file_write(handle, 'Niexp', profni)
  s = file_write(handle, 'Tiexp', profti)
  s = file_write(handle, 'Teexp', profte)
  s = file_write(handle, 'pressure_s', profp)
  s = file_write(handle, 'Nixexp', profni[0,ypeak])
  file_close, handle

end

  
