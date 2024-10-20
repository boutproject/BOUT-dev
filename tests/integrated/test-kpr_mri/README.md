test-kpr_mri
===========

 Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:

    [u]' = [ G  e ] [(-1+u^2-r)/(2u)] + [      r'(t)/(2u)        ]
    [v]    [ e -1 ] [(-2+v^2-s)/(2v)]   [ s'(t)/(2*sqrt(2+s(t))) ]
         = [ fs(t,u,v) ]
           [ ff(t,u,v) ]

 where r(t) = 0.5*cos(t),  s(t) = cos(w*t),  0 < t < 5.

 This problem has analytical solution given by
    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).

 We use the parameters:
   e = 0.5 (fast/slow coupling strength) [default]
   G = -100 (stiffness at slow time scale) [default]
   w = 100  (time-scale separation factor) [default]

 The stiffness of the slow time scale is essentially determined
 by G, for |G| > 50 it is 'stiff' and ideally suited to a
 multirate method that is implicit at the slow time scale.
