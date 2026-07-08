#!/bin/bash
set -e

PROJECT_ROOT="@PROJECT_BINARY_DIR@"
export PYTHONPATH="@PROJECT_BINARY_DIR@/tools/pylib:@PROJECT_SOURCE_DIR@/tools/pylib:$PYTHONPATH"

# Pre-build the project to prevent concurrent CMake race conditions
cmake --build "$PROJECT_ROOT"

PYTEST_ARGS=("$@")
if [ "$CI" == "true" ]; then
    TEST_FLAGS="-v"
else
    TEST_FLAGS="-q"
fi

# Use the Python executable that CMake discovered during configuration
"@Python3_EXECUTABLE@" -m pytest -m "not serial" --cache-clear -n auto --dist=loadgroup $TEST_FLAGS "$PROJECT_ROOT/tests/integrated" "${PYTEST_ARGS[@]}"
"@Python3_EXECUTABLE@" -m pytest -m serial $TEST_FLAGS "$PROJECT_ROOT/tests/integrated" "${PYTEST_ARGS[@]}"
