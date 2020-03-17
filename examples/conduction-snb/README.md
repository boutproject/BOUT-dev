SNB non-local heat conduction
=============================

The Shurtz, Nicolai and Busquet (SNB) model [Physics of Plasmas 7, 4238 (2000)](https://doi.org/10.1063/1.1289512) is a
multigroup diffusion model which captures some features of kinetic electron heat
conduction. 

If the thermal mean free path is more than a few percent of the gradient scale
length, then the collisional Spitzer-Harm (Braginskii) heat conduction model
tends to over-estimate the heat flux. Flux limiters can in some cases capture
this, but need to be tuned and in general the tuning needed varies in time and
space. In other regions the collisional model can under-estimate the heat flux
due to the high energy tail of the distribution function. This "pre-heat" can be
to some extent captured by the SNB model.

Sinusoidal perturbation
-----------------------

    python sinusoid.py
    
This test case puts a small variation in temperature on top of a periodic background. The background temperature is varied, and the SNB heat flux compared against Spitzer-Harm for a range of mean free paths.


Temperature step
----------------

    python step.py
    
This case has a large change in temperature, from nearly 1keV to 200eV over 0.3mm,
relevant to a laser-plasma experiment. This is compared against a kinetic VFP solution.

