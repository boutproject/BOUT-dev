;-------------------------------------------------------------
;+
; NAME:
;       FACTOR
; PURPOSE:
;       Find prime factors of a given number.
; CATEGORY:
; CALLING SEQUENCE:
;       factor, x, p, n
; INPUTS:
;       x = Number to factor.            in
; KEYWORD PARAMETERS:
; OUTPUTS:
;       p = Array of prime numbers.      out
;       n = Count of each element of p.  out
; COMMON BLOCKS:
; NOTES:
;       Note: see also prime, numfactors, print_fact.
; MODIFICATION HISTORY:
;       R. Sterner.  4 Oct, 1988.
;       RES 25 Oct, 1990 --- converted to IDL V2.
;       Johns Hopkins University Applied Physics Laboratory.
;
; Copyright (C) 1988, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------
 
	pro factor, x, p, n, help=hlp
 
	if (n_params(0) lt 3) or keyword_set(hlp) then begin
	  print,' Find prime factors of a given number.'
	  print,' factor, x, p, n'
	  print,'   x = Number to factor.            in'
	  print,'   p = Array of prime numbers.      out'
	  print,'   n = Count of each element of p.  out'
	  print,' Note: see also prime, numfactors, print_fact.'
	  return
	endif
 
	s = sqrt(x)			; Only need primes up to sqrt(x).
	g = fix(50 + 0.13457*s)		; Upper limit of # primes up to s.
	p = prime(g)			; Find g primes.
	n = intarr(n_elements(p))	; Divisor count.
 
	t = long(x)			; Working number.
	i = 0				; Index of test prime.
 
loop:	pt = p(i)			; Pull test prime.
	t2 = long(t/pt)			; Result after division.
	if t eq t2*pt then begin	; Check if it divides.
	  n(i) = n(i) + 1		; Yes, count it.
	  t = t2			; Result after division.
	  if t2 eq 1 then goto, done	; Check if done.
	  goto, loop			; Continue.
	endif else begin
	  i = i + 1			; Try next prime.
	  if i ge g then goto, last	; Nothing up to sqrt works.
	  goto, loop			; Continue.
	endelse
 
last:	p = [p,t]			; Residue was > sqrt, must be prime.
	n = [n,1]			; Must occur only once. (else < sqrt).
 
done:	w = where(n gt 0)
	n = n(w)			; Trim excess off tables.
	p = p(w)
	return
	end

