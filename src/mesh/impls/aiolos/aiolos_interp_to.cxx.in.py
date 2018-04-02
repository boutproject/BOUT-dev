#!/bin/python3

import gen_stencils

print("""
#include "aiolosmesh.hxx"
#include <interpolation.hxx>
""")

from common import license
print(license())

# Should we generate using the raw pointers or field operators?
gen_stencils.use_field_operator = False
# Should the numbers be printed as fraction or as floating numbers?
gen_stencils.useFloat = False
# Should the fractions be casted to floats to force compile time
# evaluation?
gen_stencils.staticCastFloat = True
gen_stencils.print_interp_to_code()
