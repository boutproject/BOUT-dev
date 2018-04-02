#!/bin/python3

from gen_derivs import print_init
from common import license

print(license())


print("""
#include "aiolosmesh.hxx"
#include <output.hxx>
#include <strings.h>
#include "aiolos_init.hxx"
""")

print_init()
