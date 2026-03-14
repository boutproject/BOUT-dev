#!/bin/bash


HERE=$(cd "$(dirname "$0")" && pwd)
export PYTHONPATH="$HERE/../tools/pylib"

pytest -m "not serial" --cache-clear -n auto --dist=loadgroup integrated $@
pytest -m serial integrated $@
