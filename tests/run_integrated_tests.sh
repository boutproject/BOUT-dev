#!/bin/bash
set -e

PROJECT_ROOT="@PROJECT_BINARY_DIR@"

export PYTHONPATH="@PROJECT_BINARY_DIR@/tools/pylib:@PROJECT_SOURCE_DIR@/tools/pylib:$PYTHONPATH"

# Pre-build the project to prevent concurrent CMake race conditions
cmake --build "$PROJECT_ROOT"

# Use the Python executable that CMake discovered during configuration
"@Python3_EXECUTABLE@" -m pytest -m "not serial" --cache-clear -n auto --dist=loadgroup -q "$PROJECT_ROOT/tests/integrated" "$@"
"@Python3_EXECUTABLE@" -m pytest -m serial "$PROJECT_ROOT/tests/integrated" "$@"
