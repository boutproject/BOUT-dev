#!/bin/python3

from gen_derivs import print_init_header
from common import license

print(license())

print("""#pragma once

#include "aiolosmesh.hxx"

""")

print_init_header()
