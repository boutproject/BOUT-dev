#!/bin/python3

import gen_derivs
#from gen_derivs import generate_index_functions_stag, generate_index_functions
from common import license

print(license())

print("""
#include "aiolosmesh.hxx"
#include "aiolos_init.hxx"
#include "aiolos_stencils.hxx"

""")

gen_derivs.generate_index_functions_stag(gen_derivs.func_tables)
gen_derivs.generate_index_functions(header_only=False)
