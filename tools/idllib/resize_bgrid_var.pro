;***************************
; resize_bgrid_var.pro
;***************************
;
; Created:  
;   09/18/2012 I. Joseph
;
;******************************

;*****************************
; FUNCTION resize_bgrid_var
;*****************************
; Linear interpolation of variable tin to a new grid size: nxout, nyout
;
; Input:
;   gin = BOUT grid data structure
;   tin = input variable
;   nxout = new radial size
;   nyout = new poloidal size
;
; Output:
;   tout
;
; Created:  
;  05/09/2012 I. Joseph
;
;******************************

function resize_bgrid_var, gin, tin, nxout=nxout, nyout=nyout 
  ix = findgen(nxout)*float(gin.nx-1)/float(nxout-1)
  iy = findgen(nyout)*float(gin.ny-1)/float(nyout-1)
    dim = size(tin)
    case dim[0] of
      0: tout = tin
      1: begin
           case dim[1] of
             gin.nx: tout = interpolate(tin,ix)
             gin.ny: tout = interpolate(tin,iy)
             else: tout = tin
           endcase
         end
      2: begin
           tout = interpolate(tin,ix,iy,/grid)
         end
      3: begin
           tout = make_array(nxout,nyout,dim[3],type=dim[4])
           for iz=0,dim[3]-1 do tout[*,*,iz] = interpolate(reform(tin[*,*,iz]),ix,iy,/grid)
         end
    endcase
    return, tout
end
