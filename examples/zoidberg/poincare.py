#!/usr/bin/env python

# Demonstrates how to generate Poincare plots of magnetic fields

import zoidberg

# Size of the domain in y (periodic)
yperiod = 10.

# Define magnetic field
magnetic_field = zoidberg.field.StraightStellarator(I_coil=0.4, radius = 1.0, yperiod = yperiod)

# Make Poincare plot
zoidberg.plot.plot_poincare(magnetic_field, 
                            np.linspace(0,0.5,5),   # x starting positions
                            0.0,                    # z starting positions
                            yperiod,                # Periodicity of y domain
                            interactive=True)       # Click on plot to add points
