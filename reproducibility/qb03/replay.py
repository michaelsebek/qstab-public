#!/usr/bin/env python3
"""Replay the two article certificates in a temporary copy, using only Python 3.8+."""
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib
import json
import re
import shutil
import subprocess
import sys

def main():
    if sys.version_info < (3,8):
        raise SystemExit("Use Python 3.8 or newer.")
    base = Path(__file__).resolve().parent
    expected = json.loads((base / "sources.json").read_text(encoding="utf-8"))
    allowed = {"mechanical/verify_mechanical.py", "mechanical/mechanical_checks.json",
               "local/check_local_certificate.py", "local/verify_local_certificate.py",
               "local/local_certificate.json", "local/local_certificate_replay.json",
               "local/local_certificate_summary.json"}
    if set(expected["sha256"]) != allowed:
        raise SystemExit("Unexpected paths in the article source inventory.")
    with TemporaryDirectory(prefix="qstab-qb03-") as work:
        dst = Path(work)
        for name, digest in expected["sha256"].items():
            src = base / name
            if not re.fullmatch(r"[0-9a-f]{64}",digest) or hashlib.sha256(src.read_bytes()).hexdigest() != digest:
                raise SystemExit("Article source changed: " + name)
            target = dst / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, target)
        for name, args in [
            ("mechanical/verify_mechanical.py", []),
            ("local/check_local_certificate.py", [str(dst/"local/local_certificate.json")])]:
            subprocess.run([sys.executable, "-I", str(dst/name), *args],cwd=dst, check=True)
    print("Both exact article checks completed. The distributed files were not modified.")
if __name__ == "__main__":
    main()
