#!/usr/bin/env python3
"""Exact rational, whole-half-line certificates for qB01 radius brackets.

Consumes the exact same rational coefficient data as benchmark_qb01.py.
Numerical references only propose bracket endpoints and a rational frequency;
every reported bracket is accepted by exact polynomial tests over Q.
No sampling or floating-point root decision is used by the certifier.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
x, w, s = sp.symbols('x w s', real=True)

def require(condition, message):
    if not bool(condition):
        raise ValueError(message)

def mul(a, b):
    a0,a1,a2,a3=a; b0,b1,b2,b3=b
    return [a0*b0-a1*b1-a2*b2-a3*b3,
            a0*b1+a1*b0+a2*b3-a3*b2,
            a0*b2-a1*b3+a2*b0+a3*b1,
            a0*b3+a1*b2-a2*b1+a3*b0]

def conj(a): return [a[0],-a[1],-a[2],-a[3]]
def norm2(a): return sp.expand(sum(z*z for z in a))
def dot(a,b): return sp.expand(sum(z*y for z,y in zip(a,b)))

def even_poly(expr):
    poly=sp.Poly(sp.expand(expr),w,domain=sp.QQ)
    require(all(k[0]%2==0 for k,c in poly.terms()), 'Frequency invariant is not even')
    return sp.Poly(sum(c*x**(k[0]//2) for k,c in poly.terms()),x,domain=sp.QQ)

def invariants(case, side):
    tau=sp.Rational(case['tau'])
    require(tau>0, 'Frequency scale must be positive')
    require(len(case['p0'])>=2 and len(case['p1'])<len(case['p0']), 'Strict degree gap in s required')
    require(any(sp.Rational(z)!=0 for z in case['p0'][-1]), 'Leading nominal coefficient is zero')
    require(any(sp.Rational(z)!=0 for a in case['p1'] for z in a), 'This dataset certifier requires a nonzero perturbation polynomial')
    require(side in ('L','R'), 'Unknown multiplication side')
    coeffs={name:[[sp.Rational(c)*tau**k for c in a]
                   for k,a in enumerate(case[name])] for name in ['p0','p1']}
    common=max(abs(c) for values in coeffs.values() for a in values for c in a)
    coeffs={name:[[c/common for c in a] for a in vals] for name,vals in coeffs.items()}
    parts=[]; full=[]
    for name in ['p0','p1']:
        vals=coeffs[name]
        xi=[sum(a[j]*(-w*w)**(k//2) for k,a in enumerate(vals) if k%2==0) for j in range(4)]
        eta=[w*sum(a[j]*(-w*w)**(k//2) for k,a in enumerate(vals) if k%2==1) for j in range(4)]
        cc=(mul(conj(eta),xi) if side=='L' else mul(xi,conj(eta)))[1:]
        parts.append((norm2(xi)+norm2(eta),cc))
        full.append([sum(a[j]*s**k for k,a in enumerate(vals)) for j in range(4)])
    (n,c0),(d,c1)=parts
    D0=even_poly(n*n-4*norm2(c0))
    D1=even_poly(d*d-4*norm2(c1))
    D01=even_poly(n*d-4*dot(c0,c1))
    # Exact independent nominal Routh-Hurwitz verification on N=sum components^2.
    nominal=sp.Poly(norm2(full[0]),s,domain=sp.QQ)
    deg=nominal.degree()
    get=lambda k: nominal.nth(k) if 0<=k<=deg else 0
    hurwitz=sp.Matrix(deg,deg,lambda i,j:get(deg-2*(j+1)+(i+1)))
    minors=[sp.factor(hurwitz[:k,:k].det()) for k in range(1,deg+1)]
    require(nominal.LC()>0 and all(z>0 for z in minors), 'Exact nominal Hurwitz check failed')
    return even_poly(n),even_poly(d),D0,D1,D01,common,minors

def primitive_positive(poly):
    content,primitive=poly.primitive()
    if content<0: primitive=-primitive
    return primitive

def main():
    inputs=json.loads((HERE/'benchmark'/'benchmark_inputs.json').read_text())
    results=json.loads((HERE/'benchmark'/'benchmark_results.json').read_text())
    refs={(r['case'],r['side']):r for r in results['references']}
    reports=[]
    for case in inputs:
        for side in case['sides']:
            n,d,D0,D1,D01,common,minors=invariants(case,side)
            ref=refs[(case['id'],side)]
            guess=sp.Rational(ref['radius80'])
            scale=sp.Integer(10)**12
            lo=sp.floor(guess*scale)/scale
            hi=(sp.floor(guess*scale)+1)/scale
            # A deliberately open extra last decimal handles an exactly rational
            # radius without presenting equality as strict robustness.
            if lo==guess: lo-=1/scale
            if hi==guess: hi+=1/scale
            Qlo=D0-2*lo**2*D01+lo**4*D1
            Tlo=n-lo**2*d
            Qprim=primitive_positive(Qlo)
            positive_roots=int(Qprim.count_roots(0,sp.oo))
            require(Tlo.eval(0)>0 and Qlo.eval(0)>0 and positive_roots==0,
                    'Exact lower-radius certificate failed for '+case['id']+'-'+side)
            tau=sp.Rational(case['tau'])
            wguess=sp.Rational(ref['omega80'])/tau
            # Frequency selected numerically; final inequality is entirely exact.
            qx=sp.floor(wguess**2*10**24)/sp.Integer(10)**24
            Qhi=D0-2*hi**2*D01+hi**4*D1
            Thi=n-hi**2*d
            det_sign=sp.sign(Qhi.eval(qx)); trace_sign=sp.sign(Thi.eval(qx))
            require(qx>=0 and (det_sign<=0 or trace_sign<=0),
                    'Exact upper-radius certificate failed for '+case['id']+'-'+side)
            row=dict(case=case['id'],side=side,
                     radius_lower=str(lo),radius_upper=str(hi),width=str(hi-lo),
                     lower_decimal=str(sp.N(lo,16)),upper_decimal=str(sp.N(hi,16)),
                     normalized_squared_frequency_witness=str(qx),
                     lower_determinant_coefficients=[str(c) for c in Qprim.all_coeffs()],
                     lower_nonnegative_root_count=positive_roots,
                     lower_trace_at_zero=str(Tlo.eval(0)),lower_det_at_zero=str(Qlo.eval(0)),
                     upper_determinant_sign=int(det_sign),upper_half_trace_sign=int(trace_sign),
                     frequency_scale=str(tau),common_polynomial_scale=str(common),
                     nominal_hurwitz_minors=[str(z) for z in minors],status='EXACT CERTIFIED')
            reports.append(row)
            print(case['id'],side,'[',row['lower_decimal'],',',row['upper_decimal'],']',row['status'])
    out={'method':'Exact rational Sturm root count for det(P0*P0-r^2 P1*P1), or the right-order Gram pair, on x=(omega/tau)^2 >=0. Positive definiteness anchored at x=0. Upper endpoint rejected at one exact rational x.',
         'statement':'Every reported interval contains the global all-direction radius. The decimal reference only proposes endpoints; exact inequalities determine acceptance.',
         'sympy':sp.__version__,'count':len(reports),'certificates':reports}
    (HERE/'global_certificates.json').write_text(json.dumps(out,indent=2)+'\n')
    keys=['case','side','lower_decimal','upper_decimal','width','lower_nonnegative_root_count','upper_determinant_sign','upper_half_trace_sign','status']
    with (HERE/'global_certificates.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys); writer.writeheader()
        writer.writerows({k:r[k] for k in keys} for r in reports)

if __name__=='__main__': main()
