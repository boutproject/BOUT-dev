# Test vector operators

DISCLAIMER: THIS IS A WORK IN PROGRESS FOLDER, EVERYTHING MAY NOT WORK AS
EXPECTED

Contains implementations used to investigate advection of vector quantities.
These quantites can pop up in for example the vorticity in dirft reduced fluid
models, where one need to find the divergence of the material derivative of the
density times the ion polarization drift.

* [1-initialCheckContravariant](/1-initialCheckContravariant/) contains a
  BOUT++ implementation which shows how to set up calculations using a vector
  field (using contravariant components), and how to plot these.
* [2-initialCheckCovariant](/2-initialCheckCovariant/) contains a BOUT++
  implementation which shows how to set up calculations using a vector field
  (using covariant components), and how to plot these.
* [3-MESChristoffel](/3-MESChristoffel/) contains a BOUT++ implementation which
  uses the Method of Exact Solution to verify the Christoffel symbols
* [4-Vpar_Grad_par](/4-Vpar_Grad_par/) shows the shows the advection of a
  vector quantity using the ```Vpar_Grad_par``` operator. NOT IMPLEMENTED
* [5-b0xGrad_dot_Grad_advection](/5-b0xGrad_dot_Grad_advection/) shows the
  advection of a vector quantity using the ```b0xGrad_dot_Grad_advection```
  operator. NOT IMPLEMENTED
* [6-V_dot_grad_advections](/6-V_dot_grad_advections/) shows the advection of a
  vector quantity using the ```V_dot_grad``` operator. NOT IMPLEMENTED

# TODOs (DELME WHEN DONE)
1. MES the Christoffel symbols
    * This will be an indicator that the implementation in mesh.cxx is correct.
    * Could be done by specifying the contravariant metric tensor, analytically
      invert it, and from that analytically calculate the Christoffel symbols
      (see 3-MESChristoffel/calculations/clebschSystemToValidate.ipynb)
        * Could be that clebschSystemToValidate.ipynb is too hard to get a
          grasp on, and needs restructure
    * Update the README in MESChristoffel when done.

2. Validate the advection operators above
    * This can be tricky, as for example MMS would only tell us that the
      derivatives are correctly implemented.
    * An idea could be visually trace a gaussian perturbation (see
      [1-initialCheckContravariant](/1-initialCheckContravariant/))
      in a magnetic field which does not spread or compress the perturbation.
      It could (on the other hand) turn out to be tricky to find a coordinate
      sytem which exercises all the Christoffel symbols.
