#!/usr/bin/env python3
"""Compare all ROCTX/NVTX region averages in rocprofv3 or nsys summaries.

Inputs may be rocprof directories (searched recursively), nsys stats text files,
nsys nvtx_sum CSV files, or directories containing nsys text reports. MPI ranks are combined with
a call-weighted mean: sum(TotalDurationNs) / sum(Calls). Speedup is the first
input's mean divided by the compared input's mean. Regions are sorted by the
first input's average duration, longest first; regions missing from that input
appear last. Total call counts can be shown for each profile with --calls. Nested ranges overlap.
"""

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path


def add(totals, name, duration, calls, path):
    duration, calls = float(duration), int(calls)
    if not math.isfinite(duration) or duration < 0 or calls < 0:
        raise ValueError(f"{path}: non-finite/negative duration or negative call count")
    totals[name][0] += duration
    totals[name][1] += calls


def region_name(name):
    # nsys prefixes ranges with their domain; the empty domain is just ':'.
    # Preserve named domains so distinct ranges do not accidentally merge.
    return name.removeprefix(":")


def read_csv(path, totals):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        columns = set(reader.fieldnames or [])
        if {"Name", "Calls", "TotalDurationNs"}.issubset(columns):
            keys = ("Name", "TotalDurationNs", "Calls")
        elif {"Range", "Instances", "Total Time (ns)"}.issubset(columns):
            keys = ("Range", "Total Time (ns)", "Instances")
        else:
            raise ValueError(
                f"{path}: expected rocprofv3 marker statistics "
                "(Name, Calls, TotalDurationNs) or nsys nvtx_sum CSV "
                "(Range, Instances, Total Time (ns))"
            )
        for row in reader:
            name, duration, calls = (row[key] for key in keys)
            add(
                totals,
                region_name(name),
                duration.replace(",", ""),
                calls.replace(",", ""),
                path,
            )


def read_text(path, totals):
    # nsys stats nvtx_sum: percentage, total ns, instances, average, median,
    # min, max, standard deviation, then the range name (possibly prefixed).
    in_table = False
    report = None
    with path.open() as stream:
        for line in stream:
            section = re.search(r"Executing '([^']+)' stats report", line)
            if section:
                report = section.group(1)
                in_table = False
                continue
            if "Total Time (ns)" in line and "Instances" in line:
                in_table = (report is None or report == "nvtx_sum") and "Range" in line
                continue
            if not line.strip():
                in_table = False
            if not in_table:
                continue
            fields = line.split(None, 8)
            if len(fields) != 9 or not re.fullmatch(r"[\d.,]+", fields[0]):
                continue
            style_and_name = fields[8].split(None, 1)
            if len(style_and_name) == 2 and style_and_name[0] in {
                "PushPop",
                "StartEnd",
            }:
                add(
                    totals,
                    region_name(style_and_name[1].strip()),
                    fields[1].replace(",", ""),
                    fields[2].replace(",", ""),
                    path,
                )


def extract(path):
    totals = defaultdict(lambda: [0.0, 0])
    if path.is_file():
        files = [path]
    elif path.is_dir():
        files = sorted(path.rglob("*marker_api_stats.csv"))
        if not files:
            files = sorted(path.rglob("*nvtx_sum*.csv"))
        if not files:
            files = sorted(
                p
                for p in path.rglob("*")
                if p.is_file()
                and (p.suffix.lower() in {".txt", ".log"} or p.name.startswith("HW3D-"))
            )
    else:
        raise ValueError(f"input does not exist: {path}")
    if not files:
        raise ValueError(
            f"{path}: no marker statistics, nvtx_sum CSV, or text reports found"
        )
    for file in files:
        if file.suffix.lower() == ".csv":
            read_csv(file, totals)
        else:
            read_text(file, totals)
    if not any(calls for _, calls in totals.values()):
        raise ValueError(
            f"{path}: no NVTX/ROCTX region statistics found (use nsys nvtx_sum or rocprofv3 marker statistics)"
        )
    return totals


def make_table(paths, profiles, labels=None, show_calls=False):
    labels = labels or [path.stem for path in paths]
    headers = ["Region"]
    for i, path in enumerate(paths):
        label = labels[i]
        if show_calls:
            headers.append(f"{label} calls")
        headers.append(f"{label} average (ms)")
        if i:
            headers.append(f"{label} speedup vs {labels[0]}")
    names = set().union(*(profile.keys() for profile in profiles))

    def sort_key(name):
        duration, calls = profiles[0].get(name, (0, 0))
        return (not calls, -duration / calls if calls else 0, name)

    ordered = sorted(names, key=sort_key)
    rows = []
    for name in ordered:
        row = [name]
        base_duration, base_calls = profiles[0].get(name, (0, 0))
        base = base_duration / base_calls if base_calls else None
        for i, profile in enumerate(profiles):
            duration, calls = profile.get(name, (0, 0))
            mean = duration / calls if calls else None
            if show_calls:
                row.append(str(calls) if calls else "N/A")
            row.append(f"{mean / 1e6:.6f}" if mean is not None else "N/A")
            if i:
                row.append(
                    f"{base / mean:.3f}x" if base is not None and mean else "N/A"
                )
        rows.append(row)
    return headers, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "profiles",
        nargs="+",
        type=Path,
        help="run directories, rocprofv3 marker CSV, or nsys nvtx_sum text/CSV; first is baseline",
    )
    parser.add_argument(
        "--csv", action="store_true", help="emit CSV instead of Markdown"
    )
    parser.add_argument(
        "--calls",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="include total call counts per profile (default: disabled)",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        metavar="LABEL",
        help="dataset labels in input order; one per profile",
    )
    args = parser.parse_args()
    if args.labels is not None and len(args.labels) != len(args.profiles):
        parser.error("--labels requires exactly one label per profile")
    labels = args.labels or [path.stem for path in args.profiles]
    try:
        profiles = [extract(path) for path in args.profiles]
        headers, rows = make_table(args.profiles, profiles, labels, args.calls)
    except (OSError, ValueError, csv.Error) as error:
        parser.error(str(error))
    if args.csv:
        writer = csv.writer(sys.stdout)
        writer.writerow(headers)
        writer.writerows(rows)
    else:
        table = [
            [
                str(cell).replace("|", "\\|").replace("\n", " ").replace("\r", " ")
                for cell in row
            ]
            for row in [headers] + rows
        ]
        widths = [max(len(row[i]) for row in table) for i in range(len(headers))]

        def format_row(row):
            return (
                "| "
                + " | ".join(cell.ljust(width) for cell, width in zip(row, widths))
                + " |"
            )

        print(format_row(table[0]))
        print("| " + " | ".join("-" * max(3, width) for width in widths) + " |")
        for row in table[1:]:
            print(format_row(row))


if __name__ == "__main__":
    main()
