"""Tests for the redocking benchmark harness."""

import json
import tempfile
import unittest

from snaclex import benchmark, provenance
from tests.fixtures import atom, component, structure


def _receptor_with_ligand():
    coords = [
        ("C", 3, 0, 0), ("O", -3, 0, 0), ("N", 0, 3, 0), ("C", 0, -3, 0),
        ("C", 2, 2, 1), ("O", -2, 2, -1), ("N", 2, -2, 1), ("C", -2, -2, -1),
        ("C", 4, 1, 0), ("C", -4, -1, 0), ("O", 1, 4, 0), ("N", -1, -4, 0),
        ("C", 3, -1, 2), ("C", -3, 1, -2), ("O", 0, 0, 4), ("N", 0, 0, -4),
    ]
    prot = [atom(el, x, y, z, name=el, res_name="LEU", res_seq=i)
            for i, (el, x, y, z) in enumerate(coords)]
    lig = component("LIG", [
        atom("C", 0.0, 0.0, 0.0, hetero=True, res_name="LIG", chain="X", res_seq=900),
        atom("C", 1.2, 0.0, 0.0, hetero=True, res_name="LIG", chain="X", res_seq=900),
        atom("O", 1.8, 1.0, 0.0, hetero=True, res_name="LIG", chain="X", res_seq=900),
        atom("N", -1.0, 0.6, 0.4, hetero=True, res_name="LIG", chain="X", res_seq=900),
        atom("C", 0.4, -1.2, 0.8, hetero=True, res_name="LIG", chain="X", res_seq=900),
        atom("O", -0.6, -0.4, -1.2, hetero=True, res_name="LIG", chain="X", res_seq=900),
    ])
    return structure(prot, [lig]), lig


class TestRedockCase(unittest.TestCase):
    def test_returns_metrics(self):
        s, lig = _receptor_with_ligand()
        out = benchmark.redock_case(s, lig, seeds=8, mc_steps=4, seed=1)
        self.assertEqual(out["n_heavy_atoms"], 6)
        self.assertIsInstance(out["rmsd"], float)
        self.assertIn("score", out)

    def test_deterministic(self):
        s, lig = _receptor_with_ligand()
        a = benchmark.redock_case(s, lig, seeds=8, mc_steps=4, seed=2)
        b = benchmark.redock_case(s, lig, seeds=8, mc_steps=4, seed=2)
        self.assertEqual(a["rmsd"], b["rmsd"])
        self.assertEqual(a["score"], b["score"])

    def test_tiny_ligand_skipped(self):
        tiny = component("ZN", [
            atom("ZN", 0, 0, 0, hetero=True, res_name="ZN", chain="X", res_seq=1),
        ])
        s = structure([atom("C", 1, 0, 0, res_seq=1)], [tiny])
        out = benchmark.redock_case(s, tiny)
        self.assertIsNone(out["rmsd"])
        self.assertIn("skipped", out)


class TestBenchmarkStructure(unittest.TestCase):
    def test_filters_by_min_heavy(self):
        s, _lig = _receptor_with_ligand()
        # The ligand has 6 heavy atoms; min_heavy=8 should exclude it.
        self.assertEqual(benchmark.benchmark_structure(s, min_heavy=8), [])
        cases = benchmark.benchmark_structure(s, min_heavy=6, seeds=8, mc_steps=4)
        self.assertEqual(len(cases), 1)


class TestSummarize(unittest.TestCase):
    def test_metrics(self):
        cases = [
            {"rmsd": 1.0}, {"rmsd": 1.5}, {"rmsd": 3.0},
            {"rmsd": None, "skipped": "x"},
        ]
        s = benchmark.summarize(cases, success_A=2.0)
        self.assertEqual(s["n_cases"], 4)
        self.assertEqual(s["n_scored"], 3)
        self.assertEqual(s["top1_success_rate"], round(2 / 3, 3))
        self.assertEqual(s["median_rmsd_A"], 1.5)

    def test_empty(self):
        s = benchmark.summarize([])
        self.assertEqual(s["n_scored"], 0)
        self.assertIsNone(s["median_rmsd_A"])


class TestRun(unittest.TestCase):
    def test_run_aggregates(self):
        s, _lig = _receptor_with_ligand()
        orig = benchmark._load_source
        benchmark._load_source = lambda src: ("FAKE", s)
        try:
            res = benchmark.run(["whatever"], seeds=8, min_heavy=6)
        finally:
            benchmark._load_source = orig
        self.assertEqual(res["summary"]["n_cases"], 1)
        self.assertEqual(res["cases"][0]["structure"], "FAKE")
        self.assertIn("run_utc", res)


class TestDockingBenchmarkProvenance(unittest.TestCase):
    def test_missing_file_returns_none(self):
        self.assertIsNone(provenance.docking_benchmark(path="/no/such/file.json"))

    def test_parses_committed_results(self):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump({
                "run_utc": "2026-06-14 00:00 UTC",
                "datasets": ["PoseBusters"],
                "summary": {"n_cases": 10, "top1_success_rate": 0.4, "median_rmsd_A": 2.3},
            }, fh)
            path = fh.name
        out = provenance.docking_benchmark(path=path)
        self.assertEqual(out["last_benchmark_utc"], "2026-06-14 00:00 UTC")
        self.assertEqual(out["top1_success_rate"], 0.4)
        self.assertEqual(out["datasets"], ["PoseBusters"])


if __name__ == "__main__":
    unittest.main()
