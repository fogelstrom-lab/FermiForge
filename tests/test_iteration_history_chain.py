import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from iteration_history_chain import concatenate_history, read_history


class HistoryChainTests(unittest.TestCase):
    def test_offsets_and_restart_preserve_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            p, c, out, out2 = [Path(tmp)/n for n in ("parent", "current", "joined", "joined2")]
            p.write_text("# iteration map_rms\n1 0.1\n2 0.05\n")
            c.write_text("# iteration map_rms\n1 0.04\n2 0.03\n")
            original = p.read_bytes(), c.read_bytes()
            self.assertEqual(concatenate_history(p,c,out), (2,2))
            _, rows = read_history(out)
            self.assertEqual([r["iteration"] for r in rows], ["1","2","3","4"])
            self.assertEqual([r["history_restart"] for r in rows], ["0","0","1","0"])
            self.assertEqual(rows[2]["local_iteration"], "1")
            self.assertEqual(rows[2]["segment"], "2")
            concatenate_history(out,c,out2)
            _, rows = read_history(out2)
            self.assertEqual(rows[4]["segment"], "3")
            self.assertEqual(rows[4]["history_restart"], "1")
            self.assertEqual(original, (p.read_bytes(), c.read_bytes()))

    def test_reject_invalid_rows_and_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            p, c, out = [Path(tmp)/n for n in ("parent", "current", "joined")]
            for content in ("# iteration a\n1 nan\n", "# iteration a\n2 1\n",
                            "# iteration a\n1\n"):
                p.write_text(content)
                with self.assertRaises(ValueError):
                    read_history(p)
            p.write_text("# iteration a\n1 1\n")
            c.write_text("# iteration b\n1 1\n")
            with self.assertRaises(ValueError):
                concatenate_history(p,c,out)


if __name__ == "__main__":
    unittest.main()
