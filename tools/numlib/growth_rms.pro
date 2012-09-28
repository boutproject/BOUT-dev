; growth_rms.pro
;
;

function deriv2, time=time, f
  if keyword_set(time) then begin
     ans = deriv(time,f)
  endif else begin
     ans = deriv(f)
  endelse  
  return, ans
end

; **************
; PRO growth_rms
; **************
;
; Calculate growth rate of data
; for average and rms values over z-coordinate
;
; Input: 
;           f        = data array: 1d-4d
;           time     = optional time array 
;   boolean average  = perform an average over x and y 
;   boolean xaverage = perform an average over x
;   boolean yaverage = perform an average over y
;
; Output:
;   rmsgrowth = growth rate of rms component (over z) 
;   dcgrowth  = growth rate of dc component (over z) 
;
; Created:  2011/08/09  I. Joseph
;
; Modified:
;
;   2012/03/28 I. Joseph: Added header
;                         Added multi-dimensional averaging capabillity
;

pro growth_rms, f, time=t, rmsgrowth=rmsgrowth, dcgrowth=dcgrowth, average=average, xaverage=xaverage, yaverage=yaverage 

  d = size(f)
  dt = d[0]
  nt = d[dt]

  if dt lt 2 then begin
    print, 'ERROR: growth_rate requires data that is 2d-4d'
    return
  endif

  avg = keyword_set(average) or keyword_set(xaverage) or keyword_set(yaverage)  
  avg_x = keyword_set(average) or keyword_set(xaverage)
  avg_y = keyword_set(average) or keyword_set(yaverage)


  if dt gt 2 and avg then begin
    case dt of
      3: begin
           favg = mean(f,dim=1)
           growth_rms, favg, time=t, rmsgrowth=rmsgrowth, dcgrowth=dcgrowth
         end
      4: begin
           if avg_y then favg = mean(f,dim=2)
           if avg_x then favg = mean(favg,dim=1)
	   growth_rms, favg, time=t, rmsgrowth=rmsgrowth, dcgrowth=dcgrowth
         end
        default: print, 'ERROR: growth_rate requires data that is 2d-4d'
    endcase  
  endif else begin

  acdc, f, rms=rms, dc=dc
  if not keyword_set(time) then time = findgen(nt)

  case dt of 
    2: begin
         rmsgrowth = deriv2(time=time, alog(rms)) 
         dcgrowth  = deriv2(time=time, alog(abs(dc))) 
       end
    3: begin
         rmsgrowth = fltarr([d[1],d[3]])
         dcgrowth  = fltarr([d[1],d[3]])
	 help, time, rmsgrowth, rms
         for i=0,d[1]-1 do begin
           rmsgrowth[i,*] = deriv2(time=time, alog(reform(rms[i,*]))) 
           dcgrowth[i,*]  = deriv2(time=time, alog(abs(reform(dc[i,*])))) 
         endfor
       end
    4: begin
         rmsgrowth = fltarr([d[1],d[2],d[4]])
         dcgrowth  = fltarr([d[1],d[2],d[4]])
	 help, time, rmsgrowth, rms
         for i=0,d[1]-1 do begin 
         for j=0,d[2]-1 do begin 
           rmsgrowth[i,j,*] = deriv2(time=time, alog(reform(rms[i,j,*]))) 
           dcgrowth[i,j,*]  = deriv2(time=time, alog(abs(reform(dc[i,j,*])))) 
         endfor
         endfor
       end
    endcase
  endelse
end
