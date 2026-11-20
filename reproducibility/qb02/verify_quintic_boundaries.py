#!/usr/bin/env python3
"""Exact continuous boundary certificates for the four IJRC interval quintics.

Requires Python >= 3.10 and SymPy >= 1.12; no MATLAB, floating arithmetic,
frequency grid, parameter grid, or numerical root solver is used.

Default: replay the retained rational Bernstein subdivision certificate.
  python verify_quintic_boundaries.py
Rebuild the certificate independently:
  python verify_quintic_boundaries.py --build

For each of six closed square faces of each family, the real degree-ten
norm polynomial is formed exactly. Its leading and constant coefficients
and its first nine principal Hurwitz determinants are certified positive
on the WHOLE face by rational tensor Bernstein coefficients. The tenth
determinant is the constant coefficient times the ninth. Thus the real
norm is strictly Hurwitz everywhere on every boundary face. This is a
finite sufficient certificate for these specific boxes, not a claim to
implement the general exact interval algorithm from the article.

The retained dyadic subdivision is only a polynomial-positivity proof
partition. It does not replace a face by finitely many sampled members.
All rational constants, including epsilon=10^-12, are retained exactly.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from math import comb
from pathlib import Path
import sys
import time

try:
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix
except ImportError:
    raise SystemExit("This exact verifier requires SymPy: python -m pip install sympy")

HERE = Path(__file__).resolve().parent
CERTIFICATE = HERE / "quintic_boundary_certificate.json"
S, X, Y = sp.symbols("s x y")
QQXY = sp.QQ.poly_ring(X, Y)
R = sp.Rational
F = Fraction


def qfrac(q):
    return F(int(q.numerator), int(q.denominator))


def base_components(kind):
    """Exact fixed quartic components, before the small real shift."""
    if kind == "O":
        return [S**4+R(1,10)*S**3+R(1,4)*S**2-R(1,10)*S-R(1,8),
                R(5,4)*S**3+R(1,5)*S**2+R(5,8)*S,
                R(3,40)*(S**2+1), sp.Integer(0)]
    return [S**4+R(3,2)*S**3-R(1,4)*S**2-R(1,2)*S-1,
            R(17,8)*S**3+2*S**2+R(15,8)*S,
            -R(3,2)*(S**3+S), sp.Integer(0)]


def components(family, t):
    """Quaternion components for normalized parameters t_j in [-1,1]."""
    kind = "O" if family in ("OOO", "EEE") else "E"
    delta = R(1,4000000) if kind == "O" else R(1,125000000)
    eps, radius = R(1,10**12), R(1,1000)
    p = [sp.expand(c.subs(S, S-delta)) for c in base_components(kind)]
    if kind == "O":
        for j in range(3):
            p[j+1] += radius*t[j]*S
    else:
        p[1] += radius*t[0]
        p[2] += radius*t[1]*S
        p[3] += radius*t[2]*S
    if family in ("OOO", "EOO"):
        p[0] += eps*S**5
    else:
        p = [sp.expand(S*c) for c in p]
        p[0] -= eps
    return p


def polynomial_terms(p):
    """Canonical exact rational coefficients for a polynomial in x,y."""
    return [(tuple(int(k) for k in mon), qfrac(c)) for mon,c in p.terms()]


def polynomial_hash(terms):
    data = [[list(e),str(c.numerator),str(c.denominator)] for e,c in terms]
    return hashlib.sha256(json.dumps(data,separators=(",",":"),ensure_ascii=True).encode()).hexdigest()


def bernstein(terms):
    """Tensor Bernstein form on [-1,1]^2, with exact rational arithmetic."""
    m = max((e[0] for e,c in terms), default=0)
    n = max((e[1] for e,c in terms), default=0)
    power = [[F(0) for _ in range(n+1)] for _ in range(m+1)]
    for (p,q),a in terms:
        for i in range(p+1):
            ax = a*comb(p,i)*(-1)**(p-i)*2**i
            for j in range(q+1):
                power[i][j] += ax*comb(q,j)*(-1)**(q-j)*2**j
    left = [[sum(power[i][j]*F(comb(k,i),comb(m,i)) for i in range(k+1))
             for j in range(n+1)] for k in range(m+1)]
    return [[sum(left[k][j]*F(comb(l,j),comb(n,j)) for j in range(l+1))
             for l in range(n+1)] for k in range(m+1)]


def split_line(line):
    left, right = [line[0]], [line[-1]]
    while len(line)>1:
        line = [(a+b)/2 for a,b in zip(line[:-1],line[1:])]
        left.append(line[0]); right.append(line[-1])
    return left, right[::-1]


def split_tensor(b, axis):
    if axis == 1:
        pairs = [split_line(row) for row in b]
        return [p[0] for p in pairs], [p[1] for p in pairs]
    cols = list(map(list,zip(*b)))
    lo,hi = split_tensor(cols,1)
    return list(map(list,zip(*lo))), list(map(list,zip(*hi)))


def minimum(b):
    return min(c for row in b for c in row)


def certify(initial, max_depth=24, max_cells=250000):
    """Generate an exhaustive disjoint dyadic leaf partition."""
    todo = [("",initial)]; leaves=[]; nodes=0
    while todo:
        path,b = todo.pop(); nodes += 1
        lower=minimum(b)
        if lower > 0:
            leaves.append({"path":path,"lower":str(lower)})
            continue
        if len(path)>=max_depth or nodes>=max_cells:
            raise RuntimeError("Bernstein budget exhausted; no positivity conclusion")
        # Alternate x,y; constant directions need not be split.
        axis=len(path)%2
        if (len(b)==1 and axis==0) or (len(b[0])==1 and axis==1): axis=1-axis
        lo,hi=split_tensor(b,axis)
        # Each symbol records the coordinate and half explicitly.
        todo.append((path+("R" if axis==0 else "U"),hi))
        todo.append((path+("L" if axis==0 else "D"),lo))
    return leaves,nodes


def replay(initial, leaves):
    """Verify positivity AND complete face coverage, not merely listed cells."""
    tree={}
    for leaf in leaves:
        path=leaf["path"]
        if any(c not in "LRDU" for c in path): raise AssertionError("Bad certificate path")
        node=tree
        for c in path:
            if "lower" in node: raise AssertionError("Overlapping certificate leaves")
            node=node.setdefault(c,{})
        if node: raise AssertionError("Duplicate/overlapping certificate leaf")
        node["lower"]=leaf["lower"]
    todo=[(tree,initial)]; nodes=0
    while todo:
        node,b=todo.pop(); nodes+=1
        if set(node)=={"lower"}:
            low=minimum(b)
            if not low>0: raise AssertionError("Nonpositive Bernstein coefficient")
            if low!=F(node["lower"]): raise AssertionError("Retained lower bound differs")
        else:
            keys=set(node)
            if keys=={"L","R"}: axis,letters=0,("L","R")
            elif keys=={"D","U"}: axis,letters=1,("D","U")
            else: raise AssertionError("Certificate does not cover the entire face")
            lo,hi=split_tensor(b,axis)
            todo.extend([(node[letters[0]],lo),(node[letters[1]],hi)])
    return nodes


def face_polynomials(family, fixed, sign):
    free=iter((X,Y)); t=[sp.Integer(sign) if j==fixed else next(free) for j in range(3)]
    p=components(family,t)
    norm=sp.Poly(sum(c*c for c in p),S)
    if norm.degree()!=10: raise AssertionError("Degree changed")
    a=norm.all_coeffs()
    yield "leading",QQXY.from_sympy(a[0])
    yield "constant",QQXY.from_sympy(a[-1])
    for k in range(1,10):
        matrix=sp.Matrix(k,k,lambda i,j:a[2*j-i+1] if 0<=2*j-i+1<=10 else 0)
        yield "Delta%d"%k,DomainMatrix.from_Matrix(matrix).convert_to(QQXY).det()


def fixture_hash():
    t=sp.symbols("t1 t2 t3")
    rows=[]
    for fam in ("OOO","EOO","EEE","EEO"):
        rows.append([fam,[[[list(e), str(c.p), str(c.q)] for e,c in sp.Poly(p,S,*t).terms()] for p in components(fam,t)]])
    return hashlib.sha256(json.dumps(rows,separators=(",",":"),ensure_ascii=True).encode()).hexdigest()


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build",action="store_true",help="regenerate retained exact certificate")
    ap.add_argument("--certificate",type=Path,default=CERTIFICATE)
    args=ap.parse_args(argv)
    started=time.monotonic(); build=args.build
    if build:
        data={"schema":"qstab.interval-quintic-boundaries.1","arithmetic":"exact rational QQ[x,y] and tensor Bernstein","fixture_sha256":fixture_hash(),"criterion":"positive leading and constant coefficients, Delta1 through Delta9; Delta10=constant*Delta9","faces":[]}
    else:
        data=json.loads(args.certificate.read_text(encoding="utf-8"))
        if data.get("schema")!="qstab.interval-quintic-boundaries.1": raise AssertionError("Wrong certificate schema")
        if data["fixture_sha256"]!=fixture_hash(): raise AssertionError("Fixture differs from retained certificate")
    count=0; total_leaves=0; max_depth=0
    expected_faces=[(fam,fix,sign) for fam in ("OOO","EOO","EEE","EEO") for fix in range(3) for sign in (-1,1)]
    if not build and len(data["faces"])!=len(expected_faces): raise AssertionError("Wrong face count")
    for idx,(family,fixed,sign) in enumerate(expected_faces):
        label="%s:t%d=%+d/1000"%(family,fixed+1,sign)
        if build: face={"family":family,"fixed_coordinate":fixed,"sign":sign,"checks":[]}
        else:
            face=data["faces"][idx]
            if (face["family"],face["fixed_coordinate"],face["sign"])!=(family,fixed,sign): raise AssertionError("Wrong face")
            if len(face["checks"])!=11: raise AssertionError("Missing Hurwitz checks")
        fleaves=0
        for j,(name,poly) in enumerate(face_polynomials(family,fixed,sign)):
            terms=polynomial_terms(poly); b=bernstein(terms)
            if build:
                leaves,nodes=certify(b)
                check={"name":name,"polynomial_sha256":polynomial_hash(terms),"bidegree":[len(b)-1,len(b[0])-1],"leaves":leaves}
                face["checks"].append(check)
            else:
                check=face["checks"][j]
                if check["name"]!=name or check["polynomial_sha256"]!=polynomial_hash(terms): raise AssertionError("Hurwitz determinant differs")
                leaves=check["leaves"]; nodes=replay(b,leaves)
            count+=1; fleaves+=len(leaves);total_leaves+=len(leaves)
            max_depth=max(max_depth,max(len(c["path"]) for c in leaves))
        if build: data["faces"].append(face)
        print("PASS "+label+" | 11 exact polynomial checks | %d proof leaves"%fleaves,flush=True)
    summary={"ok":True,"families":4,"complete_closed_faces":24,"polynomial_checks":count,"proof_leaves":total_leaves,"max_dyadic_depth":max_depth,"method":"exact rational Hurwitz determinants with full tensor-Bernstein face coverage"}
    if build:
        data["summary"]=summary
        args.certificate.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    elif data["summary"]!=summary: raise AssertionError("Summary differs")
    print(json.dumps(summary,indent=2));print("Elapsed %.3f s"%(time.monotonic()-started))
    return summary

if __name__=="__main__":
    main()
