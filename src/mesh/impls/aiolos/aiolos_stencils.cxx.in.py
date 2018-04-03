#!/bin/python3

from common import license

print("""
#include "aiolosmesh.hxx"
#include <output.hxx>
#include <interpolation.hxx>
""")

print(license())

import stencils
# Should we generate using the raw pointers or field operators?
# Does currently not work with field operator
stencils.use_field_operator = False
stencils.print_stencil_implementations()
