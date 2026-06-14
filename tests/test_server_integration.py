"""End-to-end tests against a live in-process server.

Boots the real ThreadingHTTPServer on an ephemeral port and exercises endpoints
that need no upstream network (version, static, error paths) plus the security
header pipeline.
"""

import http.client
import json
import threading
import time
import unittest
from http.server import ThreadingHTTPServer
from unittest import mock

import server


class TestServerIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.httpd = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.port = cls.httpd.server_address[1]
        cls.thread = threading.Thread(target=cls.httpd.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()

    def setUp(self):
        # Give each test a fresh rate-limit budget (shared module singletons).
        server._IP_LIMITER._buckets.clear()
        server._HEAVY_LIMITER._buckets.clear()

    def _get(self, path, headers=None):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        conn.request("GET", path, headers=headers or {})
        resp = conn.getresponse()
        body = resp.read()
        conn.close()
        return resp, body

    def _post(self, path, payload):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        body = json.dumps(payload)
        conn.request("POST", path, body=body,
                     headers={"Content-Type": "application/json"})
        resp = conn.getresponse()
        data = resp.read()
        conn.close()
        return resp, data

    def test_version_endpoint(self):
        resp, body = self._get("/api/version")
        self.assertEqual(resp.status, 200)
        data = json.loads(body)
        self.assertEqual(data["version"], server.SNACLEX_VERSION)
        self.assertTrue(data["research_only"])

    def test_security_headers_on_static(self):
        resp, _ = self._get("/")
        self.assertEqual(resp.status, 200)
        self.assertEqual(resp.getheader("X-Content-Type-Options"), "nosniff")
        self.assertEqual(resp.getheader("X-Frame-Options"), "DENY")
        self.assertIn("3Dmol.org", resp.getheader("Content-Security-Policy"))

    def test_hsts_only_behind_https_proxy(self):
        resp, _ = self._get("/api/version")
        self.assertIsNone(resp.getheader("Strict-Transport-Security"))
        resp2, _ = self._get("/api/version", {"X-Forwarded-Proto": "https"})
        self.assertIsNotNone(resp2.getheader("Strict-Transport-Security"))

    def test_static_scope_blocks_non_web_files(self):
        # server.py lives at the repo root, not under web/ — must not be served.
        resp, _ = self._get("/server.py")
        self.assertEqual(resp.status, 404)

    def test_missing_param_returns_400_json(self):
        resp, body = self._get("/api/analyze")
        self.assertEqual(resp.status, 400)
        self.assertIn("error", json.loads(body))

    # ---- async job queue (stubbed runner, no upstream network) -------
    def test_job_lifecycle_done(self):
        with mock.patch.dict(server._JOB_RUNNERS,
                             {"echo": lambda p: {"got": p}}, clear=False):
            resp, body = self._post("/api/jobs",
                                    {"kind": "echo", "params": {"hi": 1}})
            self.assertEqual(resp.status, 202)
            job_id = json.loads(body)["job_id"]

            for _ in range(50):
                r, b = self._get(f"/api/jobs/{job_id}")
                st = json.loads(b)
                if st["status"] == "done":
                    self.assertEqual(st["result"], {"got": {"hi": 1}})
                    break
                time.sleep(0.02)
            else:
                self.fail("job never completed")

    def test_job_error_is_reported(self):
        def boom(_p):
            raise ValueError("nope")

        with mock.patch.dict(server._JOB_RUNNERS, {"boom": boom}, clear=False):
            resp, body = self._post("/api/jobs", {"kind": "boom", "params": {}})
            job_id = json.loads(body)["job_id"]
            for _ in range(50):
                _r, b = self._get(f"/api/jobs/{job_id}")
                st = json.loads(b)
                if st["status"] == "error":
                    self.assertIn("nope", st["error"])
                    break
                time.sleep(0.02)
            else:
                self.fail("job error never surfaced")

    def test_unknown_job_kind_400(self):
        resp, body = self._post("/api/jobs", {"kind": "nonsense", "params": {}})
        self.assertEqual(resp.status, 400)
        self.assertIn("error", json.loads(body))

    def test_unknown_job_id_404(self):
        resp, _ = self._get("/api/jobs/deadbeef")
        self.assertEqual(resp.status, 404)

    # ---- API docs ----------------------------------------------------
    def test_api_docs(self):
        resp, body = self._get("/api/docs")
        self.assertEqual(resp.status, 200)
        data = json.loads(body)
        self.assertIn("endpoints", data)
        self.assertTrue(any(e["path"] == "/api/upload" for e in data["endpoints"]))

    # ---- structure upload --------------------------------------------
    def _post_raw(self, path, text, ctype="text/plain"):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        conn.request("POST", path, body=text, headers={"Content-Type": ctype})
        resp = conn.getresponse()
        data = resp.read()
        conn.close()
        return resp, data

    def test_upload_and_analyze(self):
        resp, body = self._post_raw("/api/upload", _SAMPLE_PDB)
        self.assertEqual(resp.status, 200)
        up = json.loads(body)
        self.assertTrue(up["upload_id"].startswith("UL"))
        self.assertGreaterEqual(up["protein_atom_count"], 10)
        # The uploaded structure is now loadable by its id like any PDB.
        resp2, body2 = self._get(f"/api/analyze?pdb={up['upload_id']}")
        self.assertEqual(resp2.status, 200)
        self.assertEqual(json.loads(body2)["protein_atom_count"],
                         up["protein_atom_count"])

    def test_upload_rejects_junk(self):
        resp, body = self._post_raw("/api/upload", "not a structure at all\n")
        self.assertEqual(resp.status, 400)
        self.assertIn("error", json.loads(body))

    def test_upload_accepts_mmcif(self):
        from tests.fixtures import mmcif_text
        rows = [{"element": "C", "name": f"C{i}", "seq": i, "x": float(i)}
                for i in range(12)]
        body_text = mmcif_text(rows)
        resp, body = self._post_raw("/api/upload", body_text)
        self.assertEqual(resp.status, 200)
        up = json.loads(body)
        self.assertTrue(up["upload_id"].startswith("UL"))
        # The viewer payload is re-serialized to PDB format, not raw mmCIF.
        self.assertNotIn("_atom_site.", up["pdb_data"])
        self.assertIn("ATOM", up["pdb_data"])

    def test_upload_rejects_oversized(self):
        from tests.fixtures import mmcif_text
        rows = [{"element": "C", "name": f"C{i}", "seq": i} for i in range(12)]
        with mock.patch.object(server, "MAX_STRUCTURE_ATOMS", 5):
            resp, body = self._post_raw("/api/upload", mmcif_text(rows))
        self.assertEqual(resp.status, 413)
        self.assertIn("too large", json.loads(body)["error"])

    # ---- large assembly -> per-chain loading -------------------------
    def _two_chain_structure(self):
        from tests.fixtures import atom, structure
        prot = (
            [atom("C", float(i), 0, 0, name=f"A{i}", res_name="LEU", chain="A", res_seq=i)
             for i in range(6)]
            + [atom("C", float(i), 5, 0, name=f"B{i}", res_name="LEU", chain="B", res_seq=i)
               for i in range(6)]
        )
        return structure(prot)

    def test_analyze_too_large_offers_chains(self):
        s = self._two_chain_structure()
        meta = {"pdb_id": "TST1", "title": "Test assembly"}
        with mock.patch.object(server, "_load_full", return_value=("ATOMS", s, meta)), \
             mock.patch.object(server, "MAX_STRUCTURE_ATOMS", 8):
            resp, body = self._get("/api/analyze?pdb=TST1")
        data = json.loads(body)
        self.assertEqual(resp.status, 200)
        self.assertTrue(data["too_large"])
        self.assertEqual(data["n_atoms"], 12)
        self.assertEqual({c["chain"] for c in data["chains"]}, {"A", "B"})

    def test_analyze_single_chain_loads(self):
        s = self._two_chain_structure()
        meta = {"pdb_id": "TST2", "title": "Test assembly"}
        with mock.patch.object(server, "_load_full", return_value=("ATOMS", s, meta)), \
             mock.patch.object(server, "MAX_STRUCTURE_ATOMS", 8):
            resp, body = self._get("/api/analyze?pdb=TST2&chain=A")
        data = json.loads(body)
        self.assertEqual(resp.status, 200)
        self.assertEqual(data["id"], "TST2-A")          # synthetic subset id
        self.assertEqual(data["protein_atom_count"], 6)  # only chain A
        self.assertNotIn("_atom_site.", data["pdb_data"])
        # The subset is now loadable by its id like any structure.
        resp2, _ = self._get("/api/analyze?pdb=TST2-A")
        self.assertEqual(resp2.status, 200)

    def test_analyze_unknown_chain_404(self):
        s = self._two_chain_structure()
        meta = {"pdb_id": "TST3", "title": "x"}
        with mock.patch.object(server, "_load_full", return_value=("ATOMS", s, meta)):
            resp, _ = self._get("/api/analyze?pdb=TST3&chain=Z")
        self.assertEqual(resp.status, 404)


def _pdb_line(rec, serial, name, res, chain, seq, x, y, z, el):
    return (f"{rec:<6}{serial:>5} {name:<4} {res:>3} {chain}{seq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}")


# A minimal but column-correct PDB with 12 protein atoms.
_SAMPLE_PDB = "\n".join(
    _pdb_line("ATOM", i + 1, el, "LEU", "A", i + 1, float(i), 0.0, 0.0, el)
    for i, el in enumerate(["N", "C", "C", "O", "C", "C",
                            "N", "C", "C", "O", "C", "C"])
) + "\nEND\n"


if __name__ == "__main__":
    unittest.main()
