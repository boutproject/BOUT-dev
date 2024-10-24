test-kpr_mri
===========

 Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:

    [f]' = [ G  e ] [(-1+f^2-r)/(2f)] + [      r'(t)/(2f)        ]
    [g]    [ e -1 ] [(-2+g^2-s)/(2g)]   [ s'(t)/(2*sqrt(2+s(t))) ]
         = [ fs(t,f,g) ]
           [ ff(t,f,g) ]

 where r(t) = 0.5*cos(t),  s(t) = cos(w*t),  0 < t < 5.

 This problem has analytical solution given by
    f(t) = sqrt(1+r(t)),  g(t) = sqrt(2+s(t)).

 We use the parameters:
   e = 0.5 (fast/slow coupling strength) [default]
   G = -100 (stiffness at slow time scale) [default]
   w = 100  (time-scale separation factor) [default]

 The stiffness of the slow time scale is essentially determined
 by G, for |G| > 50 it is 'stiff' and ideally suited to a
 multirate method that is implicit at the slow time scale.

MRI implementations of the functions are as follows:

The slow explicit RHS function:
     [-0.5 * sin(t)/(2 * f)]
     [          0          ]

The slow implicit RHS function:
      [G e] * [(-1 + f^2 - 0.5 * cos(t))/(2 * f) ]
      [0 0]   [(-2 + g^2 - cos(w * t))/(2 * g)   ]

The fast implicit RHS function:
     [0  0] * [(-1 + f^2 - 0.5 * cos(t))/(2 * f)] + [         0         ]
     [e -1]   [(-2 + g^2 - cos(w * t))/(2 * g)  ]   [-w * sin(w * t)/(2 * sqrt(2 + cos(w*t)))]
