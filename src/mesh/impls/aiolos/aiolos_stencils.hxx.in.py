#!/bin/python3

print("""
// To be included from aiolosmesh.hxx
""")

import stencils
stencils.print_stencil_implementations(header_only=True)
