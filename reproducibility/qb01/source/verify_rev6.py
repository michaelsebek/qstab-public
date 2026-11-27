#!/usr/bin/env python3
"""New exact checks for the physical E1 interpretation and revision-5 E5 variant.
Run from any directory with SymPy. Historical benchmark files are not modified.
"""
import json, platform, sys
from pathlib import Path
import sympy as S
from directional import exact_directional as e
s,g=e.s,e.gamma
Z=S.Integer(0); O=S.Integer(1); rat=S.Rational
p0=(s*s+2*s+1,Z,O,Z)
K=(Z,s,Z,O)
plant=(s*s+2*s+1,-s,O,-O)
D=(Z,O,Z,O)
checks=[]
def ck(name,cond):
    if not bool(cond): raise AssertionError(name)
    checks.append(name)
def equal(q,p): return all(S.simplify(a-b)==0 for a,b in zip(q,p))
def norm(p):return S.expand(sum(x*x for x in p))
def L(q):
    return S.Matrix.hstack(*[S.Matrix(e.qmul(q,tuple(S.Integer(i==j) for i in range(4)))) for j in range(4)])
def H(p):return e.hurwitz_minors(norm(p))
# Derive the actual closed loop from physical plant, controller and calibration.
ck('plant + controller equals E1 nominal polynomial',equal(e.qadd(plant,K),p0))
ck('nominal plant Hurwitz',all(x>0 for x in H(plant)))
ck('nominal closed loop Hurwitz',all(x>0 for x in H(p0)))
xx=S.symbols('x0:4',real=True); t=tuple(xx)
a=e.qadd(e.ONE,t)
ck('real calibration error = L(t)',L(a)-S.eye(4)==L(t))
ck('induced norm identity',S.simplify(L(t).T*L(t)-sum(v*v for v in t)*S.eye(4))==S.zeros(4))
for j,v in enumerate([(0,-1,0,0),(0,0,1,-1)]):
    ck('mechanical skew coupling '+str(j),L(v).T==-L(v))
for side in ['L','R']:
    q=e.qmul(a,K) if side=='L' else e.qmul(K,a)
    actual=e.qadd(plant,q)
    perturb=e.qmul(t,K) if side=='L' else e.qmul(K,t)
    intended=e.qadd(p0,perturb)
    ck(side+' hypercomplex feedback identity',equal(actual,intended))
    G=L(a)*L(K) if side=='L' else L(K)*L(a)
    ck(side+' real feedback matrix identity',S.simplify(L(plant)+G-L(intended))==S.zeros(4))
# Both exact radii are known from rev5; check the frequency functions and global minima independently.
w=S.symbols('w', nonnegative=True)
I=S.I
def chi(q):
    a,b,c,d=q
    return S.Matrix([[a+I*b,c+I*d],[-c+I*d,a-I*b]])
P0=chi(p0).subs(s,I*w);P1=chi(K).subs(s,I*w)
# The largest generalized eigenvalue is the reciprocal of rho^2.
# Verify the stated rho is a root and that it is the smaller contact-squared root.
for side,rr in [('L',2+w**4/(1+w)**2),('R',(w**4+2*w**2-4*w+2)/(1+w)**2)]:
    A0=P0.conjugate().T*P0 if side=='L' else P0*P0.conjugate().T
    A1=P1.conjugate().T*P1 if side=='L' else P1*P1.conjugate().T
    B=S.simplify(A0-rr*A1)
    ck(side+' frequency formula determinant',S.factor(B.det())==0)
    # Here P1 has singular values w+1 and |w-1|; trace of B has nonnegative closed form.
    tr=S.factor(S.trace(B));
    trace_expected=4*w**5/(w+1)**2 if side=='L' else 4*w*(w*w+2)**2/(w+1)**2
    ck(side+' frequency formula PSD trace',S.simplify(tr-trace_expected)==0 and trace_expected.is_nonnegative)
    ck(side+' positive radius expression',rr.subs(w,0)==2)
    # rr is positive for all w>=0 by nonsingularity of P0 and this PSD equality.
    if side=='L': ck('L global minimum certificate',S.simplify(rr-2-w**4/(1+w)**2)==0)
    else:
        b=S.sqrt(3)-1;rrstar=14-8*S.sqrt(3)
        ck('R exact nonnegative factor certificate',S.simplify((1+w)**2*(rr-rrstar)-(w-b)**2*(w*w+2*b*w+2*S.sqrt(3)))==0)
        ck('R minimum attained',S.simplify(rr.subs(w,b)-rrstar)==0)
# Exact ray witnesses and instability inside the selected calibration ball.
results={}
for side,first in [('L',O),('R',2-S.sqrt(3))]:
    term=e.qmul(D,K) if side=='L' else e.qmul(K,D)
    p=e.qadd(p0,tuple(g*v for v in term))
    minors=H(p);h=S.factor(minors[-2]);a1=rat(1,3)
    expected=8*(g-3)*(g-1)*(g*g-2*g+2) if side=='L' else 24*(g*g-4*g+1)*(g*g-2*g+2)
    ck(side+' guardian factorization',S.expand(h-expected)==0)
    ck(side+' first positive guardian root',S.simplify(S.Poly(h,g).real_roots()[0]-first)==0)
    pc=tuple(S.simplify(v.subs(g,first)) for v in p)
    lam=(Z,Z,Z,Z) if side=='L' else (Z,Z,-(S.sqrt(3)-1),Z)
    ck(side+' original polynomial contact',e.qzero(e.qeval(pc,lam,'right' if side=='L' else 'left')))
    vals=[S.factor(v.subs(g,a1)) for v in minors]
    ck(side+' at 1/3 correct stability classification',all(v>0 for v in vals) if side=='L' else vals[2]<0 and all(v!=0 for v in vals))
    results[side]={'ray_guardian':str(h),'first_gamma':str(first),'norm_radius_squared':str(S.simplify(2*first**2)), 'at_gamma_1_3_hurwitz_minors':list(map(str,vals)),'original_contact':list(map(str,lam))}
ck('chosen witness inside half-unit ball',rat(2,9)<rat(1,4))
ck('actuator radius above half-unit ball',2>rat(1,4))
ck('sensor radius below half-unit ball',S.simplify(14-8*S.sqrt(3)<rat(1,4)))
ck('right third minor exact',results['R']['at_gamma_1_3_hurwitz_minors'][2]=='-208/27')
# Pure common scale error (real t): the two families coincide and h is constant.
realp=e.qadd(p0,tuple(g*v for v in K))
realminors=H(realp)
ck('real gain norm exact',S.expand(norm(realp)-((s+1)**4+1+g*g*(s*s+1)))==0)
ck('real scale reduced Hurwitz determinant constant',realminors[-2]==48)
results['real_scale_error']={'guardian':str(realminors[-2]),'radius':'infinity','norm':str(norm(realp))}
# Newly added noncommuting E5 variant from author rev5: no floating-point decisions.
pert=(-s-1,s,-s-1,Z);p=e.qadd(p0,tuple(g*v for v in pert));minors=H(p);h=S.factor(minors[-2]); cubic=3*g**3-9*g**2+10*g-6
ck('E5 variant guardian',S.expand(h+8*(g-1)**2*cubic)==0)
ck('E5 variant cubic increasing',S.expand(S.diff(cubic,g)-(9*(g-1)**2+1))==0)
ck('E5 variant next root > 1',cubic.subs(g,1)<0)
ck('E5 variant next root rational bracket',cubic.subs(g,rat(1747,1000))<0<cubic.subs(g,rat(1748,1000)))
ck('E5 variant exact stable after contact',all(v.subs(g,rat(3,2))>0 for v in minors))
ck('E5 variant contact factor',equal(tuple(v.subs(g,1) for v in p),(s*s+s,s,-s,Z)))
ck('E5 variant zero on original polynomial',e.qzero(e.qeval(tuple(v.subs(g,1) for v in p),e.ZERO,'right')))
ck('E5 variant coefficients have no common slice',not equal(e.qmul((0,0,1,0),(0,1,-1,0)),e.qmul((0,1,-1,0),(0,0,1,0))))
report={'status':'PASS','scope':'new exact checks, no new floating-point benchmark or MATLAB/qstab run','python':platform.python_version(),'sympy':S.__version__,'number_checks':len(checks),'checks':checks,'E1':results,'E5_variant':{'p1_components':list(map(str,pert)),'guardian':str(h),'first_gamma':'1','multiplicity':2,'next_root_interval':['1747/1000','1748/1000'],'stable_at_gamma_3_2':True,'original_zero':'0'}}
out=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('rev6_exact_checks.json');out.write_text(json.dumps(report,indent=2)+'\n')
print('PASS',len(checks),'exact checks;',out)
