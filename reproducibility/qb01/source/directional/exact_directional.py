#!/usr/bin/env python3
"""Exact rational ray verification for the qB01 polynomial article.

Inputs are four real-polynomial components of p0,p1 and a nonzero rational
quaternion d.  The ray parameter is gamma in t=gamma*d.  Thus the quaternion
norm at first loss is |d|*gamma, not gamma unless |d|=1.

The computation uses exact rational/algebraic arithmetic throughout:
1. form the original left/right polynomial and its real norm;
2. compute the reduced Hurwitz determinant h(gamma);
3. square-free factor h, isolate every real root, retain multiplicities;
4. lift each positive candidate to the ORIGINAL norm via even/odd gcd;
5. reconstruct and substitute actual left/right evaluation zero witnesses.

This is a reproducible small-case symbolic verifier, not a scalable solver
benchmark or a replacement for the public qstab implementation.
"""

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import sympy as S

s, gamma, x = S.symbols("s gamma x", real=True)
ZERO = (S.Integer(0),) * 4
ONE = (S.Integer(1), S.Integer(0), S.Integer(0), S.Integer(0))
CHECKS = []


def check(name, condition):
    passed = bool(condition)
    CHECKS.append({"name": name, "passed": passed})
    if not passed:
        raise AssertionError(name)


def qadd(a, b):
    return tuple(S.expand(v + w) for v, w in zip(a, b))


def qneg(a):
    return tuple(-v for v in a)


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return tuple(S.expand(v) for v in (
        a0*b0-a1*b1-a2*b2-a3*b3,
        a0*b1+a1*b0+a2*b3-a3*b2,
        a0*b2-a1*b3+a2*b0+a3*b1,
        a0*b3+a1*b2-a2*b1+a3*b0,
    ))


def qconj(a):
    return (a[0], -a[1], -a[2], -a[3])


def qnorm2(a):
    return S.factor(sum(v*v for v in a))


def qinv(a):
    den = qnorm2(a)
    check("nonzero quaternion inversion", den != 0)
    return tuple(S.simplify(v / den) for v in qconj(a))


def qzero(a):
    return all(S.simplify(v) == 0 for v in a)


def qeval(p, q, convention):
    """p_r(q)=sum a_k q^k, p_l(q)=sum q^k a_k."""
    n = max(S.degree(v, s) if v != 0 else 0 for v in p)
    out, power = ZERO, ONE
    for k in range(int(n) + 1):
        coeff = tuple(S.expand(v).coeff(s, k) for v in p)
        term = qmul(coeff, power) if convention == "right" else qmul(power, coeff)
        out = qadd(out, term)
        power = qmul(power, q)
    return tuple(S.simplify(v) for v in out)


def hurwitz(g, size=None):
    n = int(S.degree(g, s))
    coeff = [S.expand(g).coeff(s, k) for k in range(n + 1)]
    size = n - 1 if size is None else size
    return S.Matrix(size, size, lambda i, j:
                    coeff[n - 2*(j+1) + (i+1)]
                    if 0 <= n - 2*(j+1) + (i+1) <= n else 0)


def hurwitz_minors(g):
    n = int(S.degree(g, s))
    H = hurwitz(g, n)
    return [S.factor(H[:k, :k].det()) for k in range(1, n+1)]


def positive_exact(v):
    v = S.simplify(v)
    if v.is_positive is not None:
        return v.is_positive
    answer = S.ask(S.Q.positive(v))
    if answer is not None:
        return answer
    raise ValueError(f"Exact positivity undecided: {v}")


def stringify(v):
    return S.sstr(S.factor(v))


def even_odd(g):
    n = int(S.degree(g, s))
    c = [S.expand(g).coeff(s, k) for k in range(n + 1)]
    E = S.expand(sum((-1)**k*c[2*k]*x**k for k in range(n//2+1)))
    O = S.expand(sum((-1)**k*c[2*k+1]*x**k for k in range((n+1)//2)))
    return E, O


def lift_candidate(p, N, a, tag):
    """No numerical root is used in candidate acceptance."""
    pa = tuple(S.simplify(v.subs(gamma, a)) for v in p)
    Na = S.factor(N.subs(gamma, a))
    check(tag + ": original norm identity after specialization",
          S.simplify(qnorm2(pa) - Na) == 0)
    E, O = even_odd(Na)
    Eg = S.Poly(E, x, extension=True)
    Og = S.Poly(O, x, extension=True)
    gcd = S.polys.polytools.gcd(Eg, Og).monic()
    witnesses = []

    if S.simplify(Na.subs(s, 0)) == 0:
        check(tag + ": zero frequency original polynomial", qzero(tuple(v.subs(s, 0) for v in pa)))
        witnesses.append({"omega_squared": "0", "kind": "zero frequency",
                          "right_evaluation_zero": ["0"]*4,
                          "left_evaluation_zero": ["0"]*4})

    # real_roots over exact algebraic fields retains the specified real embedding.
    for xx, multiplicity in gcd.real_roots(multiple=False):
        if not positive_exact(xx):
            continue
        check(tag + ": exact even root", S.simplify(E.subs(x, xx)) == 0)
        check(tag + ": exact odd root", S.simplify(O.subs(x, xx)) == 0)
        remainder = [S.rem(S.Poly(v, s, extension=True),
                           S.Poly(s*s+xx, s, extension=True)).as_expr()
                     for v in pa]
        B = tuple(S.simplify(S.expand(v).coeff(s, 1)) for v in remainder)
        A = tuple(S.simplify(S.expand(v).coeff(s, 0)) for v in remainder)
        check(tag + ": remainder norm relation", S.simplify(qnorm2(A)-xx*qnorm2(B)) == 0)
        check(tag + ": remainder orthogonality", S.simplify(qmul(A, qconj(B))[0]) == 0)

        if qzero(B):
            check(tag + ": spherical remainder", qzero(A))
            qr = (S.Integer(0), S.sqrt(xx), S.Integer(0), S.Integer(0))
            ql = (S.Integer(0), S.Integer(0), S.sqrt(xx), S.Integer(0))
            kind = "full imaginary sphere (remainder identically zero)"
        else:
            qr = tuple(S.simplify(v) for v in qneg(qmul(qinv(B), A)))
            ql = tuple(S.simplify(v) for v in qneg(qmul(A, qinv(B))))
            kind = "unique zero per evaluation convention on this sphere"

        check(tag + ": right witness pure imaginary", S.simplify(qr[0]) == 0)
        check(tag + ": left witness pure imaginary", S.simplify(ql[0]) == 0)
        check(tag + ": right witness norm", S.simplify(qnorm2(qr)-xx) == 0)
        check(tag + ": left witness norm", S.simplify(qnorm2(ql)-xx) == 0)
        check(tag + ": original RIGHT evaluation", qzero(qeval(pa, qr, "right")))
        check(tag + ": original LEFT evaluation", qzero(qeval(pa, ql, "left")))
        witnesses.append({
            "omega_squared": stringify(xx),
            "gcd_multiplicity": multiplicity,
            "kind": kind,
            "remainder_B": list(map(stringify, B)),
            "remainder_A": list(map(stringify, A)),
            "right_evaluation_zero": list(map(stringify, qr)),
            "left_evaluation_zero": list(map(stringify, ql)),
        })

    return {"gamma": stringify(a), "original_norm": stringify(Na),
            "even_polynomial": stringify(E), "odd_polynomial": stringify(O),
            "even_odd_gcd": stringify(gcd.as_expr()),
            "imaginary_contact": bool(witnesses), "witnesses": witnesses}


def run_case(name, p0, p1, direction, side, expected_gamma, expected_first_multiplicity=None):
    started = time.perf_counter()
    d = tuple(map(S.sympify, direction))
    p0 = tuple(map(S.sympify, p0))
    p1 = tuple(map(S.sympify, p1))
    check(name + ": rational direction", all(v.is_Rational for v in d))
    check(name + ": nonzero direction", qnorm2(d) != 0)
    check(name + ": strict degree gap", max(S.degree(v,s) if v != 0 else -1 for v in p1)
          < max(S.degree(v,s) if v != 0 else -1 for v in p0))
    perturb = qmul(d, p1) if side == "L" else qmul(p1, d)
    p = qadd(p0, tuple(gamma*v for v in perturb))
    N = S.expand(qnorm2(p))
    product = qmul(p, qconj(p))
    check(name + ": product is real norm", S.simplify(product[0]-N)==0 and qzero(product[1:]+(0,)))
    degree = int(S.degree(N, s))
    check(name + ": constant positive leading coefficient",
          not S.expand(N).coeff(s,degree).has(gamma) and positive_exact(S.expand(N).coeff(s,degree)))
    minors = hurwitz_minors(N)
    check(name + ": exact nominal Hurwitz", all(positive_exact(v.subs(gamma,0)) for v in minors))
    h = S.factor(minors[-2])
    check(name + ": reduced determinant match", S.simplify(h-hurwitz(N).det()) == 0)
    check(name + ": nonzero guardian", h != 0 and positive_exact(h.subs(gamma,0)))
    hp = S.Poly(h,gamma,domain=S.QQ)
    sf = hp.sqf_part()
    intervals = hp.intervals(eps=S.Rational(1,10)**24)
    real_roots = hp.real_roots(multiple=False)
    positive = [(a,m) for a,m in real_roots if positive_exact(a)]
    check(name + ": multiplicities preserved", sum(m for _,m in real_roots) == sum(m for _,m in intervals))
    check(name + ": square-free locations preserved", len(sf.real_roots()) == len(real_roots))
    check(name + ": candidate cardinality bound", len(positive) <= hp.degree())
    contacts = []
    for k, (a,multiplicity) in enumerate(positive):
        contact = lift_candidate(p,N,a,name+f" candidate {k+1}")
        contact["guardian_multiplicity"] = multiplicity
        contacts.append(contact)

    if positive:
        first, first_mult = positive[0]
        check(name + ": expected first root", S.simplify(first-expected_gamma)==0)
        check(name + ": first root has original contact", contacts[0]["imaginary_contact"])
        check(name + ": exact stable point before first", all(positive_exact(v.subs(gamma,first/2)) for v in minors))
        if expected_first_multiplicity is not None:
            check(name + ": expected multiplicity", first_mult == expected_first_multiplicity)
        sample_after = (first+positive[1][0])/2 if len(positive)>1 else first+1
        stable_after = all(positive_exact(v.subs(gamma,sample_after)) for v in minors)
        first_margin = stringify(S.sqrt(qnorm2(d))*first)
    else:
        check(name + ": expected infinite ray margin", expected_gamma is S.oo)
        first, first_mult, first_margin = S.oo, None, "oo"
        stable_after = None

    return {
        "name": name, "side": side, "parameterization": "t=gamma*d, gamma>=0",
        "direction_d": list(map(stringify,d)), "direction_norm": stringify(S.sqrt(qnorm2(d))),
        "p0_components": list(map(stringify,p0)), "p1_components": list(map(stringify,p1)),
        "norm_polynomial": stringify(N), "norm_degree": degree,
        "guardian": stringify(h), "guardian_degree": hp.degree(),
        "square_free_guardian": stringify(sf.as_expr()),
        "square_free_factors_with_multiplicities": [[stringify(f.as_expr()),m] for f,m in hp.sqf_list()[1]],
        "real_root_isolating_intervals": [
            {"lower": stringify(interval[0]), "upper": stringify(interval[1]), "multiplicity":m}
            for interval,m in intervals],
        "first_gamma": stringify(first), "first_quaternion_norm": first_margin,
        "first_guardian_multiplicity": first_mult,
        "stable_in_interval_immediately_after_first": stable_after,
        "candidates": contacts,
        "elapsed_seconds_symbolic_run": time.perf_counter()-started,
    }


def main():
    cases = []
    input_path=Path(__file__).resolve().parent.parent/'benchmark'/'benchmark_inputs.json'
    inputs={c['id']:c for c in json.loads(input_path.read_text())}
    benchmark_rays=[
        ('S','side_separation_left','L',(0,1,0,1),S.Integer(1),None),
        ('S','side_separation_right','R',(0,1,0,1),2-S.sqrt(3),None),
        ('R','rank_one_frequency_peak','L',(-3,1,0,0),S.Rational(1,4),None),
        ('D','repeated_singular_value_rational_ray','L',(-1,-1,0,0),(S.sqrt(5)-1)/2,None),
        ('N','near_repeated_singular_value_real_ray','L',(-1,0,0,0),1-S.Rational(1,10)**20,None),
        ('F','nearly_flat_peak_exact_critical_ray','L',(-2,-S.Rational(1,5000),0,0),S.Integer(1),None),
        ('T','quaternion_zero_tangency','L',(1,0,0,0),S.Integer(1),2),
        ('K','scaled_side_separation_left','L',(0,1,0,1),S.Integer(1),None),
        ('K','scaled_side_separation_right','R',(0,1,0,1),2-S.sqrt(3),None),
    ]
    for code,name,side,d,expected,multiplicity in benchmark_rays:
        source=inputs[code]
        p0=tuple(sum(S.Rational(row[j])*s**k for k,row in enumerate(source['p0'])) for j in range(4))
        p1=tuple(sum(S.Rational(row[j])*s**k for k,row in enumerate(source['p1'])) for j in range(4))
        case=run_case(name,p0,p1,d,side,expected,multiplicity)
        case['benchmark_case']=code
        case['same_input_source']=str(input_path)
        cases.append(case)
    cases.append(run_case("cubic_nonzero_tangency",(s**3+2*s*s+2*s+3,0,0,0),
                          (s*s+s+6,0,0,0),(1,0,0,0),"L",S.Integer(1),8))
    cases.append(run_case("infinite_ray_margin",(s+1,0,0,0),(1,0,0,0),(0,1,0,0),"L",S.oo))
    cases.append(run_case("guardian_extraneous_later",(s*s+3*s+2,0,0,0),
                          (-s-1,0,0,0),(1,0,0,0),"L",S.Integer(2)))
    cases.append(run_case("spherical_contact_odd_part_zero",(s*s+2*s+1,0,0,0),
                          (-2*s,0,0,0),(1,0,0,0),"L",S.Integer(1)))
    later = next(c for c in cases if c["name"]=="guardian_extraneous_later")["candidates"]
    check("later opposite real roots are rejected", len(later)==2 and later[0]["imaginary_contact"]
          and not later[1]["imaginary_contact"] and later[1]["gamma"]=="3")
    base={c['side']:c for c in cases if c.get('benchmark_case')=='S'}
    for c in cases:
        if c.get('benchmark_case')=='K':
            original=base[c['side']]
            check(c['name']+": exact unchanged ray margin",c['first_gamma']==original['first_gamma'])
            check(c['name']+": exact guardian scaling",
                  S.simplify(S.sympify(c['guardian'],locals={'gamma':gamma})/
                             S.sympify(original['guardian'],locals={'gamma':gamma})-S.Integer(10)**444)==0)
    result={"status":"OK", "arithmetic":"exact rational/algebraic SymPy; no floating point acceptance",
            "python_version":platform.python_version(),"sympy_version":S.__version__,
            "number_cases":len(cases),"number_checks":len(CHECKS),
            "checks":CHECKS,"cases":cases}
    output=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name("exact_directional_checks.json")
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"cases":len(cases),"checks":len(CHECKS),"output":str(output)}))


if __name__=="__main__":
    main()
