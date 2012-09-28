function extmult, a, b 
  return, reform (a # b, n_elements(a)*n_elements(b))
end

function mysub, a, b

  if a le 0 or b le 0 then begin
    mysub = min(a, b)
  endif else begin
    g=gcd(a, b, p=p, m=m)
    print, 'primes ', p
    print, 'mult   ', m
    mysub = 1
    if g gt 1 then begin
      for i = 0, n_elements(p)-1 do begin
        n = indgen(m[i],1)+1
;        print, n
;        print, p[i]^n
        ans = [mysub, extmult(mysub,p[i]^n)]
        mysub = ans(sort(ans)) 
        print, mysub
      endfor
    endif
  endelse
  return, mysub
end


function nproc, nx, ny, nleg, verbose=verbose
;Currently only works for nx a power of 2 and ny a power of 2

; NX,NY       = # grid points in the X and Y directions

; NXPE,NYPE   = # processors in the X and Y directions

; MXSUB,  MYSUB   = # grid points/processor
; MXG,    MYG     = # guard cells/processor
; MXSUB+2*MXG, MYSUB+2*MYG = # total points handled by each processor
  
; The general relation is 
;   NX = NXPE * MXSUB + 2*MXG
;   NY = NYPE * MYSUB  

; Restrictions
;   The branch cuts defined by the X-point must lie on a processor boundary
;   This implies that both NLEG and NCORE must be a multiple of MYSUB
;   

;X direction
  ;if alog(nx/2,2) NE 0 then begin
  ;  print, "nx must be a power of 2"
  ;end

  ax = alog(nx)/alog(2)
  nxp = 2^indgen(ax)
  mxp = nx/nxp

if N_PARAMS() lt 3 then nleg = 0

;Y direction

  ncore = ny-2*nleg
  ;if alog(ncore/nleg, 2) NE 0 then begin
  ;  print, "nx must be a power of 2"
  ;end

  if nleg gt 1 then begin
    ay = alog(nleg)/alog(2)
    myp = reverse(2^indgen(ay)*2)
    nyp = ny/myp
  endif else begin
    myp = 2^nleg
    nyp = ny/myp
  endelse

  np = nxp # nyp
  nxpe=nxp # (1+0*nyp)
  nype=(1+0*nxp) # nyp
  mxpe=mxp # (1+0*myp)
  mype=(1+0*mxp) # myp



  nproc = {nx:nx, ny:ny, nleg:nleg, np:np,nxpe:nxpe,mxpe:mxpe,nype:nype,mype:mype}

  if keyword_set(verbose) then begin
    print, 'NP  ',  nproc.np
    print, 'NXPE', nproc.nxpe
    print, 'NYPE', nproc.nype
    print, 'MXPE', nproc.mxpe
    print, 'MYPE', nproc.mype
  endif

  return, nproc
end



function nproc2, nx, ny, nleg, verbose=verbose
; Here nx is counted without gaurd cells nx -> NX-2*MXG

  factor,nx,px,mx 
  if n_elements(px) gt 1 then begin
    nxpe = [1, px^mx]
  endif else begin
    nxpe = 1
  endelse
  mxpe = nx/nxpe

  ncore = ny-2*nleg
  nype=mysub(ncore, nleg)  
  mype = ny/nype

  np = extmult(nxpe,nype)
  
  nproc = {nx:nx, ny:ny, nleg:nleg, np:np, nxpe:nxpe, mxpe:mxpe, nype:nype, mype:mype}

  if keyword_set(verbose) then begin
    print, 'NP  ',  nproc.np
    print, 'NXPE', nproc.nxpe
    print, 'NYPE', nproc.nype
    print, 'MXPE', nproc.mxpe
    print, 'MYPE', nproc.mype
  endif

  return, nproc
end

