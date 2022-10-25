#!/usr/bin/env python3

from sys import argv

from xbout import open_boutdataset

ds = open_boutdataset(
    argv[1] + "/BOUT.dmp.*.nc", keep_yboundaries=1, keep_xboundaries=1
)

ds.bout.animate_list(["rhs", "rhs_check", "error", "f"], animate_over="y", show=True)
