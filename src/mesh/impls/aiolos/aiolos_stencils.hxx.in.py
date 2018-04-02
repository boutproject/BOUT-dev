#!/bin/python3

from common import license

print("""
#pragma once
#include "aiolosmesh.hxx"
""")

print(license())

import gen_stencils
gen_stencils.gen_functions_normal(header_only=True)
