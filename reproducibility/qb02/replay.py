#!/usr/bin/env python3
"""Replay both independent exact interval-example certificates.

Requires Python >= 3.10 and SymPy >= 1.12. Run from any working directory:
    python path/to/qstab_public/reproducibility/qb02/replay.py
Optional fresh machine-readable report (the chosen file is overwritten):
    python path/to/replay.py --output path/to/my_exact_report.json

This does not run MATLAB and does not claim a fresh native MATLAB test.
The optional MATLAB/Symbolic replay is verify_qstab_interval_certificates.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    try:
        import sympy
    except ImportError:
        raise SystemExit('Install the optional exact-arithmetic dependency: python -m pip install sympy')
    started=time.monotonic()
    manifest=json.loads((HERE/'sources.json').read_text(encoding='utf-8'))
    for file in manifest['files']:
        path=HERE/file['path']
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if digest!=file['sha256']:
            raise SystemExit('Source/certificate SHA-256 mismatch: '+file['path'])
    results=[]
    for filename in ['verify_exact_examples.py','verify_quintic_boundaries.py']:
        result=subprocess.run([sys.executable,str(HERE/filename)],cwd=str(HERE),text=True,capture_output=True)
        print(result.stdout,end='')
        if result.stderr: print(result.stderr,end='',file=sys.stderr)
        results.append({'script':filename,'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
        if result.returncode: raise SystemExit(result.returncode)
    report={'schema':'qstab.interval-exact-replay.1','ok':True,'computedThisRun':True,
            'nativeMATLABRun':False,'arithmetic':'exact rational QQ and tensor Bernstein',
            'timestampUTC':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),
            'sympy':sympy.__version__,'elapsedSeconds':time.monotonic()-started,
            'coverage':{'algebraicChecks':28,'quarticBernsteinInequalities':56,'completeQuinticFaces':24,
                        'quinticPolynomialChecks':264,'quinticProofLeaves':370,'maximumDyadicDepth':5},
            'sourceManifest':manifest,'runs':results}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Fresh exact replay report:',args.output.resolve())
    print('OVERALL_OK | 28 algebraic checks | 56 quartic inequalities | 24 complete quintic faces')
    return report

if __name__=='__main__': main()
