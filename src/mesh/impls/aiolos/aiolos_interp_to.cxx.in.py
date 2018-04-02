#!/bin/python3

import stencils

print("""
#include "aiolosmesh.hxx"
#include <interpolation.hxx>
""")

from common import license
print(license())

# Should we generate using the raw pointers or field operators?
stencils.use_field_operator = False
# Should the numbers be printed as fraction or as floating numbers?
stencils.useFloat = False
# Should the fractions be casted to floats to force compile time
# evaluation?
stencils.staticCastFloat = True
stencils.print_interp_to_code()
