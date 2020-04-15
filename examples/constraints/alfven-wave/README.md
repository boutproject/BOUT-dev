Alfvén wave
===========

This is a simple model of an Alfvén wave, evolving the vorticity and parallel
electromagnetic potential, which calculates the electrostatic potential by
constraining `Laplace_perp(phi) - voricity` to be zero:

```
    // Calculate parallel current from Apar
    jpar = Delp2(Apar / (0.5*beta_e));

    // Electrostatic potential, calculate by constraining
    //     Laplace_perp(phi) - Vort
    // to be zero
    ddt(phi) = Delp2(phi) - Vort;

    // Vorticity equation
    ddt(Vort) = Div_par(jpar);
    
    // Parallel electric field
    ddt(Apar) = Grad_par(phi);
```
