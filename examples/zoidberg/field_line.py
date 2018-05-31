#!/usr/bin/env python

# Demonstrates how to generate 3D plots of magnetic field lines

import zoidberg

# Size of the domain in y (periodic)
yperiod = 10.

# Define magnetic field
magnetic_field = zoidberg.field.StraightStellarator(I_coil=0.4, radius = 1.0, yperiod = yperiod)

# Make 3D plot
zoidberg.plot.plot_3d_field_line(magnetic_field, 
                                 0.3,                    # x starting position
                                 0.0,                    # z starting position
                                 yperiod,                # Periodicity of y domain
                                 cycles = 20)
