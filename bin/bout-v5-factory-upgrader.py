#!/usr/bin/env python3

import boutupgrader

import argparse
import copy
import re


# Dictionary of factory methods that may need updating
factories = {
    "Interpolation": {
        "factory_name": "InterpolationFactory",
        "type_name": "Interpolation",
        "create_method": "create",
    },
    "InvertPar": {
        "factory_name": "InvertPar",
        "type_name": "InvertPar",
        "create_method": "create",
        "old_create_method": "Create",
        "arguments_changed": True,
    },
    "Mesh": {"factory_name": "Mesh", "type_name": "Mesh", "create_method": "Create"},
    "Laplacian": {
        "factory_name": "Laplacian",
        "type_name": "Laplacian",
        "create_method": "create",
    },
    "LaplaceXZ": {
        "factory_name": "LaplaceXZ",
        "type_name": "LaplaceXZ",
        "create_method": "create",
    },
    "SolverFactory": {
        "factory_name": "SolverFactory",
        "type_name": "Solver",
        "create_method": "createSolver",
    },
    "Solver": {
        "factory_name": "Solver",
        "type_name": "Solver",
        "create_method": "create",
    },
}


def find_factory_calls(factory, source):
    """Find all the places where the factory creation method is called,
    and return a list of the variable names

    Parameters
    ----------
    factory
        Dictionary containing 'factory_name' and 'create_method'
    source
        Text to search

    """
    return re.findall(
        r"""
        \s*([\w_]+)                  # variable
        \s*=\s*
        {factory_name}::
        .*{create_method}.*
        """.format(
            **factory
        ),
        source,
        re.VERBOSE,
    )


def find_type_pointers(factory, source):
    return re.findall(
        r"""
        \b{type_name}\s*\*\s*   # Type name and pointer
        ([\w_]+)\s*;            # Variable name
        """.format(
            **factory
        ),
        source,
        re.VERBOSE,
    )


def fix_declarations(factory, variables, source):
    """Fix the declaration of varables in source. Returns modified source

    Replaces `Type*` with either `std::unique_ptr<Type>` for
    declarations, or with `auto` for initialisations.

    Parameters
    ----------
    factory
        Dictionary of factory information
    variables
        List of variable names
    source
        Text to update

    """

    for variable in variables:
        # Declarations
        source = re.sub(
            r"""
            (.*?)(class\s*)?         # optional "class" keyword
            \b({type_name})\s*\*\s*  # Type-pointer
            ({variable_name})\s*;    # Variable
            """.format(
                type_name=factory["type_name"], variable_name=variable
            ),
            r"\1std::unique_ptr<\3> \4{nullptr};",
            source,
            flags=re.VERBOSE,
        )

        # Declarations with initialisation from factory
        source = re.sub(
            r"""
            (.*?)(class\s*)?       # optional "class" keyword
            ({type_name})\s*\*\s*  # Type-pointer
            ({variable_name})\s*   # Variable
            =\s*                   # Assignment from factory
            ({factory_name}::.*{create_method}.*);
            """.format(
                variable_name=variable, **factory
            ),
            r"\1auto \4 = \5;",
            source,
            flags=re.VERBOSE,
        )

        # Declarations with zero initialisation
        source = re.sub(
            r"""
            (.*?)(?:class\s*)?     # optional "class" keyword
            ({type_name})\s*\*\s*  # Type-pointer
            ({variable_name})\s*   # Variable
            =\s*                   # Assignment
            (0|nullptr|NULL);
            """.format(
                variable_name=variable, **factory
            ),
            r"\1std::unique_ptr<\2> \3{nullptr};",
            source,
            flags=re.VERBOSE,
        )

    return source


def fix_deletions(variables, source):
    """Remove `delete` statements of variables. Returns modified source

    Parameters
    ----------
    variables
        List of variable names
    source
        Text to update

    """

    for variable in variables:
        source = re.sub(
            r"(.*;?)\s*(delete\s*{variable})\s*;".format(variable=variable),
            r"\1",
            source,
        )

    return source


def fix_create_method(factory, source):
    """Fix change of name of factory `create` method"""

    if "old_create_method" not in factory:
        return source
    old_create_pattern = re.compile(
        r"({factory_name})\s*::\s*{old_create_method}\b".format(**factory)
    )
    if not old_create_pattern.findall(source):
        return source

    if factory.get("arguments_changed", False):
        print(
            "**WARNING** Arguments of {factory_name}::{create_method} have changed, and your current arguments may not work."
            " Please consult the documentation for the new arguments.".format(**factory)
        )

    return re.sub(
        r"({factory_name})\s*::\s*{old_create_method}\b".format(**factory),
        r"\1::{create_method}".format(**factory),
        source,
    )


def apply_fixes(factories, source, all_declarations=False):
    """Apply the various fixes for each factory to source. Returns
    modified source

    Parameters
    ----------
    factories
        Dictionary of factory properties
    source
        Text to update
    """

    modified = source

    for factory in factories.values():
        modified = fix_create_method(factory, modified)
        variables = find_factory_calls(factory, modified)
        if all_declarations:
            variables = variables + find_type_pointers(factory, modified)
        modified = fix_declarations(factory, variables, modified)
        modified = fix_deletions(variables, modified)

    return modified


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix types of factory-created objects")
    parser = boutupgrader.default_args(parser)

    parser.add_argument(
        "--all-declarations",
        "-a",
        action="store_true",
        help="Fix all declarations of factory types, not just variables created from factories",
    )

    args = parser.parse_args()

    for filename in args.files:
        with open(filename, "r") as f:
            contents = f.read()
        original = copy.deepcopy(contents)

        modified = apply_fixes(
            factories, contents, all_declarations=args.all_declarations
        )

        boutupgrader.apply_or_display_patch(
            filename, original, modified, args.patch_only, args.quiet, args.force
        )
