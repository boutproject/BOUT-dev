;***********************
; gen_vertices.pro
;***********************
;
; Created: 3/27/12 I. Joseph
;
; Use at own risk:
; 	Not all functions below are working yet
; 
;************************    

;***************************
; FUNCTION gen_domains
;***************************
; 
; Decompose grid into separate topological domains
; This is the beginning of a more general approach than gen_domains_uedge
;
;***************************
; Created: 3/27/12 I. Joseph
;*************************** 

FUNCTION gen_domains, g
  if not tag_test(g,'NPOL') then begin
    domains = gen_domains_uedge(g)
  endif else begin 
    if not tag_test(g,'IXSEPS1') then domains = gen_domains_new(g)
  endelse

  return, domains
END


;***************************
; FUNCTION gen_domains_uedge
;***************************
; 
; Decompose grid into separate topological domains
; Assumes UEDGE style mesh data
; Only 'single null' and 'one region' are working at this point
;
;***************************
; Created: 3/27/12 I. Joseph
;*************************** 


FUNCTION gen_domains_uedge, g
;  index={x:0,y:0}
  if tag_test(g,'IXSEPS2') and g.ixseps2 lt g.nx then begin
    print, 'GEN_DOMAINS: double null'

        
	xcore=indgen(g.ixseps1)
	dy = g.jyseps2_2-g.jyseps1_1
	ycore=indgen(dy) + g.jyseps1_1
  	domain[0] = {x:xcore,y:ycore}

	dx = g.nx-g.ixseps1
	xpvt=indgen(dx)+g.ixseps1
	dy = g.ny - g.jyseps2_2 
	ypvt = [indgen(dy), indgen(dy)+g.jyseps2_2]
  	domain[1] = {x:xpvt,y:ypvt}
	
	dx = g.ixsep2-g.ixseps1
	xsol=indgen(dx)+g.ixseps1
	ysol=indgen(g.ny)
  	domain[2] = {x:xsol,y:ysol}

	dx = g.nx-g.ixseps2
	xpvt=indgen(dx)+g.ixseps1
	dy = g.ny - g.jyseps2_2 
	ypvt = [indgen(dy), indgen(dy)+g.jyseps2_2]
  	domain[1] = {x:xpvt,y:ypvt}

	xsol=indgen(g.ixseps2)
	dy = g.jyseps2_2-g.jyseps1_1
	ysol=indgen(dy) + g.jyseps1_1
  	domain[0] = {x:xsol,y:ysol}
	
	dx = g.nx-g.ixseps2
	xsol=indgen(dx)+g.ixseps1
	ysol=indgen(g.ny)
  	domain[2] = {x:xsol,y:ysol}

  endif else begin
    if fix(g.nx) gt fix(g.ixseps1) then begin
      print, 'GEN_DOMAINS: single null'
;      domain = replicate(index,3)
    domain = ptrarr(3,1,/allocate_heap)

  ; CORE
      x=indgen(g.ixseps1)
      onex = 1+0*x
      dy = g.jyseps2_2-g.jyseps1_1
      y=indgen(dy) + g.jyseps1_1 + 1
      oney = 1+0*y
      x2 = x # oney
      y2 = onex # y
      core = {x:fix(x),y:fix(y),x2:fix(x2),y2:fix(y2),xbndry:[-1,2],ybndry:[0,0]} 
      domain[0] = ptr_new(core, /no_copy, /allocate_heap)

  ; PVT
      x = indgen(g.ixseps1)
      onex = 1+0*x
      dy1 = g.jyseps1_1+1 
      dy2 = g.ny-g.jyseps2_2-1
      y = [indgen(dy1), indgen(dy2)+g.jyseps2_2+1]
      oney = 1+0*y
      x2 = x # oney
      y2 = onex # y
      pvt = {x:fix(x),y:fix(y),x2:fix(x2),y2:fix(y2),xbndry:[-1,2],ybndry:[-1,-1]} 
      domain[1] = ptr_new(pvt, /no_copy, /allocate_heap)

   ; SOL	
      dx = g.nx-g.ixseps1
      x=indgen(dx)+g.ixseps1
      onex = 1+0*x
      y=indgen(g.ny)
      oney = 1+0*y
      x2 = x # oney
      y2 = onex # y
      sol = {x:fix(x),y:fix(y),x2:fix(x2),y2:fix(y2),xbndry:[1,-1],ybndry:[-1,-1]} 
      domain[2] = ptr_new(sol, /no_copy, /allocate_heap)

  endif else begin
    print, 'GEN_DOMAINS: one region'
    domain = ptrarr(3,1,/allocate_heap)
    domain = ptrarr(1,1,/allocate_heap)
    x=indgen(g.nx)
    onex = 1+0*x
    y=indgen(g.ny)
    oney = 1+0*y
    x2 = x # oney
    y2 = onex # y
    core = {x:x,y:y,x2:x2,y2:y2,xbndry:[-1,-1],ybndry:[-1,-1]} 
    domain[0] = ptr_new(core, /no_copy, /allocate_heap)

  endelse
  endelse

  
  return, domain
END

;**********************
;FUNCTION modpos, xplus, minus, yplus, yminus
;  used by gen_vertices
;**********************

FUNCTION modpos, a, b, positive=positive
  ans = a mod b
  if keyword_set(positive) and ans lt 0 then ans = ans + b
  return, ans
END


FUNCTION xplus, d, id, x0 

  nx = size((*d[id]).x,/n_elements)
  bndry = (*d[id]).xbndry[1]

  if x0+1 lt nx then begin
    xp = (*d[id]).x[x0+1]
  endif else begin
    if bndry eq id then begin
      xm = (*d[id]).x[modpos(x0+1,nx,/pos)]
    endif else begin
      if bndry lt 0 then begin
        if x0 lt 0    then xm = (*d[id]).x[0] 
        if x0 ge nx-2 then xm = (*d[id]).x[nx-1]
      endif else begin 
        xp = (*d[bndry]).x[x0+1-nx]
      endelse
    endelse
  endelse

  RETURN, xp
END 

FUNCTION xminus, d, id, x0 

  nx = size((*d[id]).x,/n_elements)
  bndry = (*d[id]).xbndry[0]

  if x0 gt 0 and x0 le nx then begin
    xm = (*d[id]).x[x0-1]
  endif else begin
    if bndry eq id then begin
      xm = (*d[id]).x[modpos(x0-1,nx,/pos)]
    endif else begin
      if bndry lt 0 then begin
        if x0 le 1  then xm = (*d[id]).x[0] 
        if x0 ge nx then xm = (*d[id]).x[nx-1]
      endif else begin 
        nb = size((*d[bndry]).x,/dimensions)
        xm = (*d[bndry]).x[nb+x0-1]
      endelse
    endelse
  endelse

  RETURN, xm
END

FUNCTION yplus, d, id, y0 

  ny = size((*d[id]).y,/n_elements)
  bndry = (*d[id]).ybndry[1]
 
  if y0+1 lt ny and y0 gt -1 then begin
    yp = (*d[id]).y[y0+1]
  endif else begin
    if bndry eq id then begin
      yp = (*d[id]).y[modpos(y0+1,ny,/pos)]
    endif else begin
      if bndry lt 0 then begin
        if y0 lt 0    then yp = (*d[id]).y[0] 
        if y0 ge ny-2 then yp = (*d[id]).y[ny-1]
      endif else begin 
        yp = (*d[bndry]).y[y0+1-ny]
      endelse
    endelse
  endelse

  RETURN, yp
END 

FUNCTION yminus, d, id, y0 

  ny = size((*d[id]).y,/n_elements)
  bndry = (*d[id]).ybndry[0]

  if y0 gt 0 and y0 lt ny then begin
    ym = (*d[id]).y[y0-1]
  endif else begin
    if bndry eq id then begin
      ym = (*d[id]).y[modpos(y0-1,ny,/pos)]
    endif else begin
      if bndry lt 0 then begin
        if y0 le 1  then ym = (*d[id]).y[0] 
        if y0 ge ny then ym = (*d[id]).y[ny-1]
      endif else begin 
        nb = size((*d[bndry]).y,/dimensions)
        ym = (*d[bndry]).y[nb+y0-1]
      endelse
    endelse
  endelse 

  RETURN, ym
END

;***********************
; PRO gen_vertices 
;***********************
;
; Generates vertices surrounding each BOUT mesh Rxy and Zxy grid point 
;
; Input:
;   g = BOUT grid data structure
;
; Optional Input: 
;   file = write gout to NETCDF file; default filename 'bout.vertices.nc'
;   append = add output to g data structure
;
; Output:
;       Output partitioned by domain index id = 0,nd-1
;       dp = pointer to array of output structures 
;           (this was used to easily traverse data; could be changed to a structure in the future)
;           inspect data via help, /str, (*dp[id])
;           x,y,x2,y2 = BOUT X,Y grid points
;           rxy,zxy   = BOUT R,Z grid centers
;           rv, zv    = BOUT R,Z grid cells (see below)
;
;       Output mapped to Rxy, Zxy grid
;       rv = r vertices: 0=center, 1...4 = vertices, ccw 
;       zv = z vertices: 0=center, 1...4 = vertices, ccw 
;       dv = domain indices mapped to rxy,zxy cell centers
;
;***************************
; Created: 3/27/12 I. Joseph
;
; Modified:
;   9/18/12 I. Joseph
;    - changed to procedure, changed some variable names
;    - added file write option
;***************************

PRO gen_vertices, g, rv=rv, zv=zv, dv=dv, dp=dp, gout=gout, append=append, file=file, time=time
  if keyword_set(time) then time0=systime(/sec)
  rv = fltarr([g.nx,g.ny,5])
  zv = fltarr([g.nx,g.ny,5])
  dv = fltarr([g.nx,g.ny])

  dp = gen_domains_uedge(g)
  nd = size(dp,/n_elements)

  for id = 0, nd-1 do begin
    print, '  domain:', id
    nx = size((*dp[id]).x,/n_elements)
    ny = size((*dp[id]).y,/n_elements) 
    rxy = fltarr(nx,ny)
    zxy = fltarr(nx,ny)
    for ix = 0, nx-1 do begin
    for iy = 0, ny-1 do begin
       jx0 = (*dp[id]).x[ix]
       jy0 = (*dp[id]).y[iy]
       rxy[ix,iy] = g.rxy[jx0,jy0]
       zxy[ix,iy] = g.zxy[jx0,jy0]
    endfor
    endfor

    xs =  fix((*dp[id]).xbndry[0] lt 0)
    xf =  fix((*dp[id]).xbndry[1] lt 0)
    ys =  fix((*dp[id]).ybndry[0] lt 0)
    yf =  fix((*dp[id]).ybndry[1] lt 0)
    print, id, nx, ny, nx-xs-xf, ny-ys-yf
    print, xs, xf, xs+xf, ys, yf, ys+yf
    r = fltarr(nx-xs-xf,ny-ys-yf,5)
    z = fltarr(nx-xs-xf,ny-ys-yf,5)
    for ix = 0, (nx-1-xs-xf) do begin
;     if id eq 0 then print, id, ix, xs, nx+xs+xf-3 
    for iy = 0, (ny-1-ys-yf) do begin

; Find neighboring grid points
       jxm = xminus(dp,id,ix+xs) ; *dp[id]).x[ix]
       jx0 = (*dp[id]).x[ix+xs]
       jxp = xplus(dp,id,ix+xs) ; (*dp[id]).x[ix+2]
       jym = yminus(dp,id,iy+ys)   ; (*dp[id]).y[iy]
       jy0 = (*dp[id]).y[iy+ys]
       jyp = yplus(dp,id,iy+ys) ; (*dp[id]).x[ix+2]

; Interpolate intermediate vertices 
       r[ix,iy,0] = g.rxy[jx0,jy0]
       z[ix,iy,0] = g.zxy[jx0,jy0]
       r[ix,iy,1] = 0.25*(g.rxy[jx0,jy0]+g.rxy[jx0,jyp]+g.rxy[jxp,jyp]+g.rxy[jxp,jy0])
       z[ix,iy,1] = 0.25*(g.zxy[jx0,jy0]+g.zxy[jx0,jyp]+g.zxy[jxp,jyp]+g.zxy[jxp,jy0])
       r[ix,iy,2] = 0.25*(g.rxy[jx0,jy0]+g.rxy[jxp,jy0]+g.rxy[jxp,jym]+g.rxy[jx0,jym])
       z[ix,iy,2] = 0.25*(g.zxy[jx0,jy0]+g.zxy[jxp,jy0]+g.zxy[jxp,jym]+g.zxy[jx0,jym])
       r[ix,iy,3] = 0.25*(g.rxy[jx0,jy0]+g.rxy[jx0,jym]+g.rxy[jxm,jym]+g.rxy[jxm,jy0])
       z[ix,iy,3] = 0.25*(g.zxy[jx0,jy0]+g.zxy[jx0,jym]+g.zxy[jxm,jym]+g.zxy[jxm,jy0])
       r[ix,iy,4] = 0.25*(g.rxy[jx0,jy0]+g.rxy[jxm,jy0]+g.rxy[jxm,jyp]+g.rxy[jx0,jyp])
       z[ix,iy,4] = 0.25*(g.zxy[jx0,jy0]+g.zxy[jxm,jy0]+g.zxy[jxm,jyp]+g.zxy[jx0,jyp])
       rv[jx0,jy0,*] = r[ix,iy,*]
       zv[jx0,jy0,*] = z[ix,iy,*]
       dv[jx0,jy0] = id + 1 
    endfor
    endfor
    (*dp[id])=create_struct((*dp[id]),'rxy',rxy,'zxy',zxy,'r',r,'z',z)
  endfor

  gout = create_struct(g,'DV',dv,'RV',rv,'ZV',zv)
  if keyword_set(append) then g = gout
  if keyword_set(file) then begin
    if size(file,/tname) ne 'STRING' then begin
      file = 'bout.vertices_'+strtrim(g.nx,2)+'x'+strtrim(g.ny,2)+'.nc'
    endif 
    print, 'gen_vertices: writing to '+file
;   Write structure in reverse order to see usual order
    s=file_export(file, reverse_struct(gout))
  endif

  if keyword_set(time) then begin
    time1=systime(/sec)
    print, 'gen_vertices: finished in '+strtrim(time1-time0,2)+' sec'
  endif 

END


;***********************
; PRO interp_rz
;***********************
;
; Interpolate to unstructured mesh
; I doubt that this is working now
;
;***************************
; Created: 3/27/12 I. Joseph
;***************************

FUNCTION interp_rz, g, fdata, Ri, Zi, domains=d, data=data
  COMMON gen_domains_com, nd, xs, xe, ys, ye

   nd = size(g.npol)
   d  = gen_domains(g)
   data  = replicate({f:0.,r:0.,z:0.,chi:0.,tri:0.},nd)

; Determine whether the points are inside the domain

  ans = fltarr(dim)
  for id=0,nd-1 do begin
    triangulate, g.d[id], g.Zxy[id], data[id].tri
    data[id].r=g.Rxy[ d[id].x, d[id].y ]
    data[id].z=g.Zxy[ d[id].x, d[id].y ]
    data[id].f=fdata[ d[id].x, d[id].y ]
    ones = 0*data[id].r + 1
    data[id].chi=griddata(data[id].r,data[id].z,ones, /linear, triangles=data[id].tri, /grid, xout=Ri[*,0], yout=Zi[0,*], missing=0)
  endfor 
 
  dim=size(Ri)
  ans=fltarr(dim)
  for i=0,dim(0)-1 do begin
  for i=0,dim(0)-1 do begin
    id = -1
    repeat begin
	id = id + 1 
    endrep until data[id].chi GT 0  
    if id le nd-1 then begin
	c[id].chi=griddata(data[id].r,data[id].z,data[id].f, /linear, triangles=data[id].tri,/grid, xout=Ri[*,0], yout=Zi[0,*], missing=0)
    endif
  endfor 
  endfor

  RETURN, ans
END


