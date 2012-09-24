;***************************
; resize_bgrid.pro
;***************************
;
; Created:  
;   03/26/2012 I. Joseph
;
; Modified:
;   09/18/2012 I. Joseph
;
;******************************



;*****************************
; FUNCTION resize_bgrid
;*****************************
; Linear interpolation of data structure of variables defined on bout grid to a new grid size: nxout, nyout
;
; Input:
;   gin = BOUT grid data structure; default structure to be resized
;   din = optional data structure to be resized
;   nxout = new radial size
;   nyout = new poloidal size
;
; Output:
;   ans = output data structure
;
; Created:  
;  05/09/2012 I. Joseph
;
;******************************

function resize_bgrid, gin, datastr=datastr, nxout=nxout, nyout=nyout

  if not keyword_set(nxout) and keyword_set(nyout) then begin
    print, 'resize_boutgrid: input parameters nxout and nyout must be specified'
    return, 0
  end

  ans = create_struct('nxout',nxout,'nyout',nyout)
  ix = findgen(nxout)*float(gin.nx-1)/float(nxout-1)
  iy = findgen(nyout)*float(gin.ny-1)/float(nyout-1)

  din = gin
  if keyword_set(datastr) then din = datastr

  tnames = tag_names(din)
  for it = 0, n_tags(din)-1 do begin
    tin = din.(it)
    tout = resize_bgrid_var(gin, tin, nxout=nxout, nyout=nyout)
    ans = create_struct(ans,tnames[it],tout)  
   endfor

   if not keyword_set(datastr) then begin
     ans.nx = nxout 
     ans.ny = nyout
     xfactor0 = float(nxout)/float(gin.nx)
     yfactor0 = float(nyout)/float(gin.ny) 
     xfactor1 = float(nxout-1)/float(gin.nx-1)
     yfactor1 = float(nyout-1)/float(gin.ny-1) 

     ans.ixseps1 = floor(gin.ixseps1*xfactor1)
     ans.ixseps2 = nxout; *******floor(gin.ixseps2*xfactor0)
     ans.ny_inner = floor(gin.ny_inner*yfactor0)
     ans.jyseps1_1 = floor((gin.jyseps1_1+1)*yfactor1-1)
     ans.jyseps1_2 = floor(gin.jyseps1_2*yfactor0)
     ans.jyseps2_1 = floor(gin.jyseps2_1*yfactor0)
     ans.jyseps2_2 = nyout - floor((gin.ny-gin.jyseps2_2-1)*yfactor1-1)
   endif

   return, ans
end



