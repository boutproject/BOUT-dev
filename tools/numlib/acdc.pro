; acdc.pro

; *************
; FUNCTION dimz
; **************
; Guess at z-dimension
;
; Input: 
;   f = data array: 1d-4d 
;   boolean time = is time included in final dimension?
;
; Output:
;   dz = dimension of z-coordinate
;        assumes that 
;
; Created:  2012/03/28  I. Joseph
;

function dimz, f, last=last

  dz  = size(f,/n_dimensions)

  undefined = dz eq 0 or (dz le 1 and not keyword_set(last)) or dz gt 4
  if undefined then begin
    print, 'ERROR: requires data that is 1d-4d'
    return, -1
  endif

  if dz eq 4 or not keyword_set(last) then dz = dz-1
  return, dz
end

; ******
; PRO dc
; ******
; Calculate average over z-coordinate
;
; Input: 
;   f    = data array: 1d-4d
;   zarray = optional array of z-coordinates   
;
; Output:
;   dc = z-average
;
; Created:  2012/03/28  I. Joseph
;

pro dc, f, zarray=zarray, dc=dc, last=last
  d  = size(f)
  dz = dimz(f,last=last)
  nz = d[dz]

  if not keyword_set(z) then begin
    dc = mean(f,dim=dz)
  endif else begin
    zlen = z[nz-1]-z[0]
    case d[0] of
      1: dc = int_tabulated(z, f)/zlen
      2: begin
           dc = fltarr(d[1])
           if not keyword_set(time_included) then begin
             for i=0, d[1]-1 do $
               dc[i] = int_tabulated(z, reform(f[i,*]))/zlen 
           endif else begin
             for t=0, d[1]-1 do $
               dc[t] = int_tabulated(z, reform(f[*,t]))/zlen 
           endelse
         end
      3: begin
           if not keyword_set(time_included) then begin
             dc = fltarr(d[1],d[2])
             for i=0, d[1]-1 do $
             for j=0, d[2]-1 do $
               dc[i,j] = int_tabulated(z, reform(f[i,j,*]))/zlen
           endif else begin
             dc = fltarr(d[1],d[3])
             for i=0, d[1]-1 do $
             for t=0, d[3]-1 do $
               dc[i,t] = int_tabulated(z, reform(f[i,*,t]))/zlen
           endelse  
         end
      4: begin
           dc = fltarr(d[1],d[2],d[4])
           for i=0, d[1]-1 do $
           for j=0, d[2]-1 do $
           for t=0, d[4]-1 do $
           dc[i,j,t] = int_tabulated(z, reform(f[i,j,*,t]))/zlen  
         end
      default: print, 'ERROR: dc: requires an f that is 1d-4d'
    endcase
  endelse
end

; ********
; PRO acdc
; ********
; Calculate fluctuations over z-coordinate
;
; Input: 
;   f      = data array: 1d-4d
;   zarray = optional array of z-coordinates  
;
; Output:
;   dc  = z-average
;   ac  = f - dc
;   rms = mean(abs(ac)) averaged over z-coordinate
;
; Created:  2012/03/28  I. Joseph
;

pro acdc, f, zarray=zarray, dc=dc, ac=ac, rms=rms, last=last
  d  = size(f)
  dz = dimz(f,last=last)
  nz = d[dz]
  ac = f
  dc, f, zarray=zarray, dc=dc, last=last
   
  case d[0] of
      1: begin 
           ac = f - dc
         end
      2: begin 
           if keyword_set(last) then begin 
             for z=0, nz-1 do ac[*,z] = ac[*,z] - dc
           endif else begin
             for z=0, nz-1 do ac[z,*] = ac[z,*] - dc
           endelse
         end 
      3: begin 
           if keyword_set(last) then begin 
             for z=0, nz-1 do ac[*,*,z] = ac[*,*,z] - dc
           endif else begin
             for z=0, nz-1 do ac[*,z,*] = ac[*,z,*] - dc
           endelse

         end
      4: begin 
           for z=0, nz-1 do ac[*,*,z,*] = ac[*,*,z,*] - dc
         end 
      default: print, 'ERROR: acdc: requires data that is 1d-4d'
  endcase

  rms = mean(abs(ac),dim=dz)
end
