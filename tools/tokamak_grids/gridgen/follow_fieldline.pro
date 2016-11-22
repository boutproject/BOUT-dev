; given a starting and ending point, returns indices along the field line between the two points
;   start_point and end_point should each be (Ri,Zi) index coordinates

pro follow_fieldline, interp_data, start_point, end_point, ri_out=ri, zi_out=zi

  ri = start_point[0]
  zi = start_point[1]

  ; value at starting point
  fval = interpolate(interp_data.f,ri,zi)
  fcut = interp_data.f
  fcut[0:1,*] = 0.0
  contour_lines, fcut, findgen(interp_data.nx), findgen(interp_data.ny), levels=fval, path_info=info, path_xy=xy
  ; weird behaviour of info, so take out important values
  NN = n_elements(xy[0,*])
  offsets = (info.offset)
  ; ri = REFORM(xy[0,offset:(offset+NN-1)])
  ; zi = REFORM(xy[1,offset:(offset+NN-1)])

  ; loop through to find the closest points to the starting and ending point
  start_diff = fltarr(NN)
  end_diff = fltarr(NN)
  for i=0,NN-1 do begin
    start_diff[i] = sqrt( (xy[0,i] - start_point[0])^2 + (xy[1,i] - start_point[1])^2 )
    end_diff[i] = sqrt( (xy[0,i] - end_point[0])^2 + (xy[1,i] - end_point[1])^2 )
  endfor
  start_ind = where(start_diff EQ min(start_diff))
  end_ind = where(end_diff EQ min(end_diff))

  if(n_elements(offsets) EQ 1) then begin
    closed = (info.type)[0]
    NN = (info.N)[0]
  endif else begin
    for i=1,n_elements(offsets)-1 do begin
      if(start_ind LT offsets[i]) then begin
        NN = (info.N)[i-1]
        xy = xy[*,offsets[i-1]:offsets[i]-1]
        start_ind -= offsets[i-1]
        end_ind -= offsets[i-1]
        closed = (info.type)[i-1]
        break
      endif
    endfor
  endelse

  if(closed) then begin
    ; check for clockwise array - calculate delta theta between two adjacent points
    R0 = mean(xy[0,*])
    Z0 = mean(xy[1,*])
    rstart = xy[0,0] - R0
    zstart = xy[1,0] - Z0
    theta_start = atan(zstart,rstart)
    rnext = xy[0,1] - R0
    znext = xy[1,1] - Z0
    theta_next = atan(znext,rnext)
    ; reverse array if it's counter-clockwise
    if(theta_next - theta_start GT 0.0) then begin
      xy = reverse(xy,2)
      start_ind = NN - 1 - start_ind
      end_ind = NN - 1 - end_ind
    endif
    ; now that it's clockwise, makesure end_ind > start_ind, otherwise shift
    if(end_ind LT start_ind) then BEGIN
      xy = shift(xy,[0,-end_ind-1])
      start_ind -= end_ind + 1
      end_ind = NN-1
    endif
  endif

  ; plot,xy[0,*],xy[1,*],psym=4,/iso
  ; oplot,xy[0,start_ind],xy[1,start_ind],psym=2
  ; stop

  ; put together array from start to end with poitns in between that lie on field line
  for i=0,abs(end_ind[0]-start_ind[0])-2 do begin
    if start_ind GT end_ind then begin
      ri = [ri,xy[0,start_ind[0]-i-1]]
      zi = [zi,xy[1,start_ind[0]-i-1]]
    endif else begin
      ri = [ri,xy[0,start_ind[0]+i+1]]
      zi = [zi,xy[1,start_ind[0]+i+1]]
    endelse
  endfor
  ri = [ri,end_point[0]]
  zi = [zi,end_point[1]]

end
