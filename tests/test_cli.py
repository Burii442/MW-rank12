"""Small failure-path checks; mathematical verification is run by verify.py."""
import importlib.util
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("verify", ROOT / "verify.py")
verify = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)


class InputChecks(unittest.TestCase):
    def test_rejects_nonrational_data(self):
        for text in ("0.5", "1/0", "I", "__import__('os')"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                verify.rational(text)

    def test_safe_latex_parser(self):
        self.assertEqual(verify.latex_scalar(r"$\frac{-4 - 3 i}{5}$"),
                         -verify.S.Rational(4, 5)-verify.S.I*verify.S.Rational(3, 5))
        with self.assertRaises(ValueError):
            verify.latex_scalar("__import__('os')")

    def test_mismatched_paper_replaces_old_pass_report(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            paper, report = folder / "paper.tex", folder / "results.json"
            paper.write_text("% BEGIN EMBEDDED NODE TABLE\n% END EMBEDDED NODE TABLE\n")
            report.write_text('{"status": "PASS"}\n')
            process = subprocess.run([sys.executable, "-O", str(ROOT / "verify.py"),
                                      "--paper", str(paper), "--output", str(report)],
                                     capture_output=True, text=True)
            self.assertEqual(process.returncode, 1, process.stdout + process.stderr)
            self.assertEqual(json.loads(report.read_text())["status"], "FAIL")
            self.assertIn("differs", process.stderr)

    def test_optimized_python_keeps_mathematical_checks(self):
        process = subprocess.run([sys.executable, "-O", "-c",
                                  "import verify; verify.require(False, 'active check')"],
                                 cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(process.returncode, 0)
        self.assertIn("active check", process.stderr)

    def test_missing_data_replaces_old_pass_report(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source = folder / "source"
            source.mkdir()
            program, report = source / "verify.py", folder / "results.json"
            program.write_bytes((ROOT / "verify.py").read_bytes())
            report.write_text('{"status": "PASS"}\n')
            process = subprocess.run([sys.executable, str(program), "--output", str(report)],
                                     capture_output=True, text=True)
            self.assertEqual(process.returncode, 1, process.stdout + process.stderr)
            result = json.loads(report.read_text())
            self.assertEqual(result["status"], "FAIL")
            self.assertIn("node_table.json", result["error"])
            self.assertNotIn("data_sha256", result)

    def test_output_cannot_replace_repository_files(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "source"
            source.mkdir()
            program = source / "verify.py"
            program.write_bytes((ROOT / "verify.py").read_bytes())
            for name in ("verify.py", "data/node_table.json", "README.md", "requirements.txt",
                         "tests/test_cli.py", ".github/workflows/verify.yml", "CITATION.cff"):
                with self.subTest(name=name):
                    target = source / name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    if target != program:
                        target.write_text("repository input\n")
                    before = target.read_bytes()
                    process = subprocess.run([sys.executable, str(program), "--output", str(target)],
                                             capture_output=True, text=True)
                    self.assertEqual(process.returncode, 2, process.stdout + process.stderr)
                    self.assertIn("outside the repository", process.stderr)
                    self.assertEqual(target.read_bytes(), before)

    def test_default_report_is_outside_repository_and_independent_of_cwd(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source, current = folder / "source", folder / "current"
            source.mkdir()
            current.mkdir()
            program = source / "verify.py"
            program.write_bytes((ROOT / "verify.py").read_bytes())
            process = subprocess.run([sys.executable, str(program)], cwd=current,
                                     capture_output=True, text=True)
            self.assertEqual(process.returncode, 1, process.stdout + process.stderr)
            report = folder / "MW-rank12-results" / "results.json"
            self.assertEqual(json.loads(report.read_text())["status"], "FAIL")
            self.assertFalse((source / "results.json").exists())
            self.assertFalse((current / "results.json").exists())

    def test_output_cannot_create_report_inside_repository(self):
        with patch.object(sys, "argv", ["verify.py", "--output", str(ROOT / "new-results.json")]), \
                patch.object(verify, "write_report") as write_report, \
                redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            verify.main()
        self.assertEqual(error.exception.code, 2)
        write_report.assert_not_called()

    def test_output_symlink_cannot_replace_repository_input(self):
        with tempfile.TemporaryDirectory() as folder:
            alias = Path(folder) / "input.json"
            try:
                alias.symlink_to(ROOT / "data" / "node_table.json")
            except OSError as error:
                self.skipTest(f"symlinks unavailable: {error}")
            with patch.object(sys, "argv", ["verify.py", "--output", str(alias)]), \
                    patch.object(verify, "write_report") as write_report, \
                    redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                verify.main()
            self.assertEqual(error.exception.code, 2)
            write_report.assert_not_called()

    def test_report_path_inside_repository_cannot_escape_through_symlink(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            source, destination = folder / "source", folder / "reports"
            source.mkdir()
            destination.mkdir()
            alias = source / "reports"
            try:
                alias.symlink_to(destination, target_is_directory=True)
            except OSError as error:
                self.skipTest(f"symlinks unavailable: {error}")
            with patch.object(verify, "HERE", source), \
                    patch.object(sys, "argv", ["verify.py", "--output", str(alias / "results.json")]), \
                    patch.object(verify, "write_report") as write_report, \
                    redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                verify.main()
            self.assertEqual(error.exception.code, 2)
            write_report.assert_not_called()

    def test_output_cannot_replace_external_manuscript(self):
        with tempfile.TemporaryDirectory() as folder:
            paper = Path(folder) / "paper.tex"
            paper.write_text("manuscript\n")
            process = subprocess.run([sys.executable, str(ROOT / "verify.py"),
                                      "--paper", str(paper), "--output", str(paper)],
                                     capture_output=True, text=True)
            self.assertEqual(process.returncode, 2)
            self.assertIn("must not overwrite the manuscript", process.stderr)
            self.assertEqual(paper.read_text(), "manuscript\n")

    def test_old_pass_is_replaced_before_computation(self):
        with tempfile.TemporaryDirectory() as folder:
            report = Path(folder) / "results.json"
            report.write_text('{"status": "PASS"}\n')

            def compute(rows):
                self.assertEqual(json.loads(report.read_text())["status"], "RUNNING")
                return {"status": "PASS", "identity_count": 0}

            with patch.object(sys, "argv", ["verify.py", "--output", str(report)]), \
                    patch.object(verify, "compute", side_effect=compute), redirect_stdout(io.StringIO()):
                self.assertEqual(verify.main(), 0)
            self.assertEqual(json.loads(report.read_text())["status"], "PASS")

    def test_keyboard_interrupt_replaces_old_pass(self):
        with tempfile.TemporaryDirectory() as folder:
            report = Path(folder) / "results.json"
            report.write_text('{"status": "PASS"}\n')

            def compute(rows):
                self.assertEqual(json.loads(report.read_text())["status"], "RUNNING")
                raise KeyboardInterrupt()

            with patch.object(sys, "argv", ["verify.py", "--output", str(report)]), \
                    patch.object(verify, "compute", side_effect=compute), redirect_stderr(io.StringIO()):
                self.assertEqual(verify.main(), 130)
            self.assertEqual(json.loads(report.read_text())["status"], "INTERRUPTED")

    def test_abrupt_exit_does_not_leave_old_pass(self):
        with tempfile.TemporaryDirectory() as folder:
            report = Path(folder) / "results.json"
            report.write_text('{"status": "PASS"}\n')
            script = ("import os, sys, verify; "
                      "verify.compute = lambda rows: os._exit(17); "
                      f"sys.argv = ['verify.py', '--output', {str(report)!r}]; verify.main()")
            process = subprocess.run([sys.executable, "-c", script], cwd=ROOT,
                                     capture_output=True, text=True)
            self.assertEqual(process.returncode, 17, process.stdout + process.stderr)
            self.assertEqual(json.loads(report.read_text())["status"], "RUNNING")


if __name__ == "__main__":
    unittest.main()
