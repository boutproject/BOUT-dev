;-------------------------------------------------------------
;+
; NAME:
;     PRIME
; PURPOSE:
;     Return an array with the specified number of prime numbers.
; EXPLANATATION:
;     This procedure is similar to PRIMES in the standard IDL distribution,
;     but stores results in a common block, and so is much faster 
;
; CALLING SEQUENCE:
;       p = prime(n)
; INPUTS:
;       n = desired number of primes, scalar positive integer
; OUTPUTS:
;       p = resulting array of primes, vector of positive integers
; COMMON BLOCKS:
;       prime_com
; NOTES:
;       Note: Primes that have been found in previous calls are
;         remembered and are not regenerated.
; MODIFICATION HISTORY:
;       R. Sterner  17 Oct, 1985.
;       R. Sterner,  5 Feb, 1993 --- fixed a bug that missed a few primes.
;       Converted to IDL V5          March 1999
;
; Copyright (C) 1985, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	function prime,n, help=hlp
 
	common prime_com, max, pmax
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Return an array with the specified number of prime numbers.'
	  print,' p = prime(n)'
	  print,'   n = desired number of primes.    in'
	  print,'   p = resulting array of primes.   out'
	  print,' Note: Primes that have been found in previous calls are'
	  print,'   remembered and are not regenerated.'
	  return, -1
	endif
 
	if n_elements(max) eq 0 then max = 0	; Make MAX defined.
	if n le max then return, pmax[0:n-1]	; Enough primes in memory.
	p = lonarr(n)				; Need to find primes.
	if max eq 0 then begin			; Have none now. Start with 8.
	  p[0] = [2,3,5,7,11,13,17,19]
	  if n le 8 then return, p[0:n-1]	; Need 8 or less.
	  i = 8					; Need more than 8.
	  t = 19L				; Search start value.
	endif else begin			; Start with old primes.
	  p[0] = pmax				; Move old primes into big arr.
	  i = max				; Current prime count.
	  t = p[max-1]				; Biggest prime so far.
	endelse
 
loop:	if i eq n then begin			; Have enough primes.
	  max = n				; Count.
	  pmax = p				; Array of primes.
	  return, p				; Return primes.
	endif
loop2:	t = t + 2				; Next test value, t.
	it = 1					; Start testing with 1st prime.
loop3:	pr = p[it]				; Pick next test prime.
	pr2 = pr*pr				; Square it.
	if pr2 gt t then begin			; Selected prime > sqrt(t)?
	  i = i + 1				; Yes, count
	  p[i-1] = t				; and store new prime.
	  goto, loop				; Go check if done.
	endif
	if pr2 eq t then goto, loop2		; Test number, t, was a square.
	if (t mod pr) eq 0 then goto, loop2	; Curr prime divides t.
	it = it + 1				; Check next prime.
	goto, loop3
	end
