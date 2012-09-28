; ***************
; FUNCTION gcd
; ***************
; Calculate greatest common denominator of 2 #'s
;
; Input: 
;   a,b    = 2 input #'s
;
; Output:
;   gcd = greatest common denominator
;   p = list of prime factors dividing both a and b
;   m = multiplicity of each prime factor
;
; Created:  2012/03/29  I. Joseph
;
; Modified: 2012/04/26  I. Joseph
;    Added header

function gcd, a, b, prime=p, multiplicity=m

  factor, a, pa, ma
  factor, b, pb, mb

  na = size(pa,/n_elements) 
  nb = size(pb,/n_elements) 
  
  nlong  = max([na, nb], ilong)
;  print, na, nb, nlong, ilong


  if ilong eq 0 then begin
    plong  = pa
    mlong  = ma
    pshort = pb
    mshort = mb
    nshort = nb
  endif else begin
    plong  = pb
    mlong  = mb
    pshort = pa
    mshort = ma
    nshort = na
  endelse

  ntotal = 0
  for is = 0,nshort-1 do begin
    ntotal = ntotal + total(plong eq pshort[is])
  endfor
  if ntotal eq 0 then begin
    p=1
    m=1
  endif else begin
    p = intarr(ntotal)
    m = intarr(ntotal)

    ctotal=0  
    for is = 0,nshort-1 do begin
      il = where(plong eq pshort[is])
      if il ge 0 then begin
        p[ctotal] = pshort[is]
        m[ctotal] = min([mshort[is],mlong[il]])
        ctotal = ctotal + 1
      endif
    endfor  
    if ctotal ne ntotal then print, 'Mismatch in ntotal'
  endelse

  return, product(p^m)
end
