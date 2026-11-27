#!/usr/bin/env python3
"""Explicit independent replay of the qB01 exact source certificates.

This optional command is never launched by the MATLAB toolbox. Original
source scripts run in a temporary copy; packaged historical evidence stays
unchanged. SymPy 1.14.0 is the recorded reference environment.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path.cwd() / 'qb01_exact_results')
    args = parser.parse_args()
    if not __debug__:
        raise RuntimeError('Run without -O: original source uses assertions.')
    import sympy
    here = Path(__file__).resolve().parent
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    report = {'schema': 'qstab.oneparameter.python-replay.1', 'ok': False,
              'python': platform.python_version(), 'sympy': sympy.__version__,
              'runtimeKind': 'independent Python exact verification, not MATLAB',
              'computedThisRun': True, 'checks': {}, 'files': [], 'errors': []}
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='qstab_qb01_exact_') as temporary:
        work = Path(temporary) / 'source'
        shutil.copytree(here / 'source', work)
        jobs = [('directional/exact_directional.py', 'directional/exact_directional_checks.json'),
                ('certify_global.py', 'global_certificates.json'),
                ('verify_rev6.py', 'rev6_exact_checks.json')]
        for script, result_name in jobs:
            print('Running', script, flush=True)
            log = out / (Path(script).stem + '.log')
            with log.open('w') as stream:
                run = subprocess.run([sys.executable, str(work / script)], cwd=work,
                                     stdout=stream, stderr=subprocess.STDOUT)
            if run.returncode:
                report['errors'].append(f'{script}: exit {run.returncode}; inspect {log.name}')
                continue
            result = work / result_name
            destination = out / result.name
            shutil.copy2(result, destination)
            data = json.loads(destination.read_text())
            if 'directional/' in result_name:
                report['checks']['directionalCases'] = data['number_cases']
                report['checks']['directionalChecks'] = data['number_checks']
                good = data['status'] == 'OK' and data['number_cases'] == 13 and data['number_checks'] == 376 and all(x['passed'] for x in data['checks'])
            elif result_name == 'global_certificates.json':
                report['checks']['globalCertificates'] = data['count']
                good = data['count'] == 9 and all(x['status'] == 'EXACT CERTIFIED' for x in data['certificates'])
            else:
                report['checks']['additionalIdentities'] = data['number_checks']
                good = data['status'] == 'PASS' and data['number_checks'] == 42
            if not good:
                report['errors'].append(f'{script}: incomplete or failed exact checks')
            report['files'].append({'name': destination.name,
                                    'sha256': hashlib.sha256(destination.read_bytes()).hexdigest()})
    report['elapsedSeconds'] = time.perf_counter() - started
    report['ok'] = not report['errors'] and len(report['files']) == 3
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
    return 0 if report['ok'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
