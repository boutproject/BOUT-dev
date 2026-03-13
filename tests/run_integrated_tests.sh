#!/bin/bash


export PYTHONPATH='../tools/pylib'

pytest -m "not serial" --cache-clear -n auto --dist=loadgroup integrated $@
pytest -m serial integrated $@
