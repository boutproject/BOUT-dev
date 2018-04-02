#!/bin/python3

from common import license

print("""
#include "aiolosmesh.hxx"
#include <output.hxx>
#include <interpolation.hxx>
""")

print(license())

import gen_stencils
# Should we generate using the raw pointers or field operators?
# Does currently not work with field operator
gen_stencils.use_field_operator = False
gen_stencils.gen_functions_normal()
