#!/bin/bash

HERE=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$HERE/../tools/pylib"

# Pre-build the project to prevent concurrent CMake race conditions
cmake --build "$HERE/.." # Compiles the parent directory (your build root)

pytest -m "not serial" --cache-clear -n auto --dist=loadgroup integrated $@
pytest -m serial integrated $@
