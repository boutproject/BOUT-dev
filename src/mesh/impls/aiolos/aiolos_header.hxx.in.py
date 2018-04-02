#!/bin/python3

from derivs import generate_index_functions
from common import license

print(license())

print(generate_index_functions(header_only=True))
