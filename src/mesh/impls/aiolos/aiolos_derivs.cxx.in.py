#!/bin/python3

import derivs
#from derivs import generate_index_functions_stag, generate_index_functions
from common import license

print(license())

print("""
#include "aiolosmesh.hxx"
#include "aiolos_init.hxx"
#include "aiolos_stencils.hxx"

""")

derivs.generate_index_functions_stag(derivs.func_tables)
derivs.generate_index_functions(header_only=False)
