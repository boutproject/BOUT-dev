#!/usr/bin/env python3
import configparser
import itertools
import os


def make_test_functions_include():
    """
    Make test_functions.cxx needed for test_initial
    """
    datadir = "data"
    inputfile = os.path.join(datadir, "BOUT.inp")

    # Read the input file
    config = configparser.ConfigParser()
    with open(inputfile, "r") as f:
        config.read_file(itertools.chain(['[global]'], f), source=inputfile)

    # Find the variables that have a "function" option
    varlist = [key for key, values in config.items() if 'function' in values]

    # Remove the coordinate arrays
    for coord in ["var_x", "var_y", "var_z"]:
        varlist.remove(coord)

    # Make the test case
    cxx_snippet = """
      Field3D {name};
      create_and_dump({name}, "{name}");
    """

    with open("test_functions.cxx", "w") as f:
        for var in varlist:
            f.write(cxx_snippet.format(name=var))


if __name__ == "__main__":
    make_test_functions_include()
