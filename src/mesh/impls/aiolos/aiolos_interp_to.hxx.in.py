#!/bin/python3

import stencils

print("""
// to be included from aiolosmesh.hxx

""")


stencils.print_interp_to_code(header_only=True)
