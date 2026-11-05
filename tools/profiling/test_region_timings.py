"""Regression checks for region comparison (run with python3 -m unittest discover)."""

import tempfile
import unittest
from pathlib import Path

from extract_region_timings import extract, make_table


class RegionTimingsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def test_weighted_rank_statistics_and_missing_regions(self):
        header = "Name,Calls,TotalDurationNs\n"
        self.write(
            "run/rank0_marker_api_stats.csv", header + "custom region,1,1000000\n"
        )
        self.write(
            "run/rank1_marker_api_stats.csv", header + "custom region,3,9000000\n"
        )
        first = extract(self.root / "run")
        other = self.write("other.csv", header + "custom region,2,2500000\nextra,1,0\n")
        headers, rows = make_table(
            [self.root / "run", other], [first, extract(other)], show_calls=True
        )
        self.assertIn("run average (ms)", headers)
        self.assertEqual(
            rows[0], ["custom region", "4", "2.500000", "2", "1.250000", "2.000x"]
        )
        self.assertEqual(rows[1], ["extra", "N/A", "N/A", "1", "0.000000", "N/A"])

    def test_nsys_text_preserves_spaces_and_domains_and_skips_other_reports(self):
        path = self.write(
            "nsys.txt",
            """Executing 'nvtx_sum' stats report
 Time (%)  Total Time (ns)  Instances  Avg (ns)  Med (ns)  Min (ns)  Max (ns)  StdDev (ns)  Style  Range
 50.0 2,000,000 2 1000000 1000000 1000000 1000000 0 PushPop :custom region
 50.0 2,000,000 1 2000000 2000000 2000000 2000000 0 StartEnd domain:other region

Executing 'cuda_gpu_kern_sum' stats report
 Time (%)  Total Time (ns)  Instances  Avg (ns)  Med (ns)  Min (ns)  Max (ns)  StdDev (ns)  Name
 100.0 9000000 1 9000000 9000000 9000000 9000000 0 unrelated kernel
""",
        )
        self.assertEqual(
            dict(extract(path)),
            {"custom region": [2000000, 2], "domain:other region": [2000000, 1]},
        )

    def test_nsys_csv_and_directory_discovery(self):
        self.write(
            "nsys/run_nvtx_sum.csv",
            'Range,Instances,Total Time (ns),Style\n:custom region,2,"2,000,000",PushPop\n',
        )
        self.assertEqual(
            dict(extract(self.root / "nsys")), {"custom region": [2000000, 2]}
        )

    def test_bad_or_unsupported_inputs_fail(self):
        for content in (
            "Name,Calls,TotalDurationNs\nx,1,nan\n",
            "Name,Calls,TotalDurationNs\nx,-1,10\n",
            "Name,Calls,TotalDurationNs\nx,0,0\n",
            "Name,Instances,Total Time (ns)\nkernel,1,100\n",
        ):
            with self.subTest(content=content), self.assertRaises(ValueError):
                extract(self.write("bad.csv", content))


if __name__ == "__main__":
    unittest.main()
