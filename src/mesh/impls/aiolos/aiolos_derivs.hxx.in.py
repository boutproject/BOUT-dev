#!/bin/python3

import derivs

print("""
//to be included from aiolosmesh.hxx

""")

derivs.generate_index_functions_stag(derivs.func_tables,header_only=True)
derivs.generate_index_functions(header_only=True)
