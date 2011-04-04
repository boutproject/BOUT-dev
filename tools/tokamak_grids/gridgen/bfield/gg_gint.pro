pro gint_load, f, n, xmin, xmax, ymin, ymax, rcentr, bcentr
;
; Load in memory 2D matrix f(n,n) and grid dimensions
; x, y are uniform grid points
;------------------------------------------------------;

 COMMON gint, g 
 DCT2D, f, ff, n ;-calculate discrete cosine transform
 g={f:f, ff:ff, n:n, xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax, rcentr:rcentr, bcentr:bcentr} 

end


function gint_fval, x0, y0, deriv=deriv
;
; Return value interpolated from pre-loaded arrays
;-------------------------------------------------------;
 COMMON gint, g

;-calculate fractional indices
 i0=(g.n-1)*(x0-g.xMin)/(g.xMax-g.xMin)
 j0=(g.n-1)*(y0-g.yMin)/(g.yMax-g.yMin)


 res=EvalCosPfast(g.ff, g.n, x0=i0, y0=j0)
 
      ;-res is defined as [f, df/dx, df/dy]
      df_dx=res[1]*(g.n-1)/(g.xMax-g.xMin)
      df_dy=res[2]*(g.n-1)/(g.yMax-g.yMin)


 if not keyword_set(DERIV) then DERIV=''


 CASE deriv of

     '': val=res[0]

     'x': val=df_dx

     'y': val=df_dy

 ENDCASE


 return, val
end



function get_psi, x, y, nonorm=nonorm, deriv=deriv
;
; Calculate flux and its derivatives at a given point[s] (x,y)
;-------------------------------------------------------------;


n=n_elements(x)
if (n ne n_elements(y)) then STOP, 'x, y dimensions not equal'

psi=dblarr(n)


for i=0,n-1 do begin
    psi[i]=Gint_Fval(x[i],y[i], deriv=deriv)
endfor


sz=SIZE(x)
if (sz[0] eq 2) then psi=REFORM(psi,sz[1],sz[2])


;if keyword_set(DERIV) then begin
;    psi=psi/(p.psixpt-p.psimagx)       
;endif else begin
;    psi=(psi-p.psimagx)/(p.psixpt-p.psimagx)       
;endelse


;
;
;
return, psi
end



function get_Br, r,z, nonorm=nonorm
 return, get_psi(r,z,deriv='y', nonorm=nonorm)
end


function get_Bz, r,z, nonorm=nonorm
 return, -get_psi(r,z,deriv='x', nonorm=nonorm)
end
