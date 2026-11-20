#!/usr/bin/env python3
"""Independent exact rational checks for IJRC interval-paper examples.

Requires Python >= 3.10 and SymPy >= 1.12. Default recomputes all claims
and compares exact rational results with the retained certificate.
Use --build to regenerate that certificate. Displayed decimal numbers
are diagnostics only; all pass/fail decisions use rational arithmetic.
The quartic polynomial bounds certify whole rectangles, not grid samples.
No MATLAB, NumPy, SciPy, numerical root solver or external package command
is invoked. Repeated/endpoint algebraic roots use SymPy exact root counts.
"""
import sympy as S
from fractions import Fraction as F
from math import comb
import argparse, hashlib, json, pathlib, time
s,v,d = S.symbols('s v d'); Q=S.Rational
root=pathlib.Path(__file__).resolve().parent
out={'schema':'qstab.interval-exact-examples.1','source':'IJRCN_manuscript_v7.tex','exact_checks':[],'bernstein':[]}
def ck(name, b, details=None):
 out['exact_checks'].append({'name':name,'pass':bool(b),'details':details}); print(name,bool(b),flush=True)
def norm(C): return S.Poly(sum(x*x for x in C),s)
def routh(poly):
 c=S.Poly(poly,s).all_coeffs(); n=len(c)-1; M=(n+2)//2
 rows=[list(c[::2])+[S.Integer(0)]*(M-len(c[::2])), list(c[1::2])+[S.Integer(0)]*(M-len(c[1::2]))]
 for i in range(2,n+1):
  a,b=rows[-2:]; assert b[0]!=0
  rows.append([S.cancel((b[0]*a[j+1]-a[0]*b[j+1])/b[0]) for j in range(M-1)]+[S.Integer(0)])
 return [x[0] for x in rows]
def cv(x): return F(int(S.numer(x)),int(S.denom(x)))
def bern(name,expr,lo,hi,dm,target=0,degrees=None):
 pp=S.Poly(expr,v,d); m=pp.degree(v); n=pp.degree(d)
 if degrees: assert m<=degrees[0] and n<=degrees[1]; m,n=degrees
 lo,hi,dm=cv(lo),cv(hi),cv(dm); span=hi-lo
 A=[[F(0) for _ in range(n+1)] for _ in range(m+1)]
 for (p,j),c in pp.terms():
  c=cv(c)*dm**j
  for i in range(p+1): A[i][j]+=c*comb(p,i)*lo**(p-i)*span**i
 B=[[sum(A[i][j]*F(comb(k,i),comb(m,i)) for i in range(k+1)) for j in range(n+1)] for k in range(m+1)]
 C=[[sum(B[k][j]*F(comb(l,j),comb(n,j)) for j in range(l+1)) for l in range(n+1)] for k in range(m+1)]
 low=min(z for row in C for z in row); target=cv(target)
 obj={'name':name,'interval':[str(lo),str(hi)],'shift_max':str(dm),'bidegree':[m,n],'coefficients':(m+1)*(n+1),'min_float':float(low),'claimed_strict_lower_bound':str(target),'pass':low>target,'minimum_exact':str(low),'polynomial_sha256':hashlib.sha256(json.dumps([[list(e),str(c.p),str(c.q)] for e,c in pp.terms()],separators=(',',':')).encode()).hexdigest()}
 out['bernstein'].append(obj); print('bern',name,float(low),'target',float(target),low>target,flush=True)
 return low
# base component polynomials, ascending powers conceptually
star=[s**4+Q(1,10)*s**3+Q(1,4)*s**2-Q(1,10)*s-Q(1,8),Q(5,4)*s**3+Q(1,5)*s**2+Q(5,8)*s,Q(3,40)*(s**2+1),S.Integer(0)]
qstar=[s**3+Q(1,10)*s**2+Q(1,2)*s,Q(1,4)*s**2+Q(1,10)*s+Q(1,8),Q(3,40)*s,Q(3,40)]
circ=[s**4+Q(3,2)*s**3-Q(1,4)*s**2-Q(1,2)*s-1,Q(17,8)*s**3+2*s**2+Q(15,8)*s,-Q(3,2)*(s**3+s),S.Integer(0)]
qcirc=[s**3+Q(3,2)*s**2+Q(7,8)*s,Q(9,8)*s**2+Q(1,2)*s+1,-Q(3,2)*s**2,-Q(3,2)*s]
for name,P,R in [('OOO',star,qstar),('EOO',circ,qcirc)]:
 ck(name+' norm factorization',norm(P).as_expr()==S.expand((s*s+1)*norm(R).as_expr()))
 rr=routh(norm(R)); ck(name+' cubic norm is Hurwitz',all(x>0 for x in rr),[str(x) for x in rr])
r=Q(1,1000); do=Q(1,4000000); de=Q(1,125000000); eps=Q(1,10**12)
for name,base,shift,which in [('OOO',star,do,'high'),('EOO',circ,de,'high'),('EEE',star,do,'low'),('EEO',circ,de,'low')]:
 C=[S.expand(x.subs(s,s-shift)) for x in base]
 if which=='high': C[0]+=eps*s**5
 else: C=[S.expand(s*x) for x in C]; C[0]-=eps
 rr=routh(norm(C)); j=9 if which=='high' else 7
 lo,hi=(-Q(1,500000),-Q(1,1000000)) if name in ['OOO','EEE'] else (-Q(7,100000000),-Q(3,50000000))
 ck(name+' quintic central Routh',all(x>0 for k,x in enumerate(rr) if k!=j) and lo<rr[j]<hi,{'index':j,'negative_entry_float':float(rr[j]),'bounds':[str(lo),str(hi)],'first_column_float':[float(x) for x in rr]})
# OOO exact tensor Bernstein certificate
print('OOO table',flush=True)
a=d**4-d**3/10-6*d*d*v+d*d/4+3*d*v/10+d/10+v*v-v/4-Q(1,8)
b=-4*d**3+3*d*d/10+4*d*v-d/2-v/10-Q(1,10)
c1=-5*d**3/4+d*d/5+15*d*v/4-5*d/8-v/5
c2=3*(d*d+1-v)/40; y1=15*d*d/4-2*d/5-5*v/4+Q(5,8); y2=-3*d/20
D=c1*c1+c2*c2; P=v*b*b-D
n=[-a*b*c1-D*y1,-a*b*c2-D*y2,S.Integer(0)]
for name,expr,l,h,t,deg in [
 ('OOO D-3/100',D-Q(3,100),Q(9,10),Q(11,10),0,(2,6)),
 ('OOO 41/1000-D',Q(41,1000)-D,Q(997,1000),Q(1003,1000),0,(2,6)),
 ('OOO a2-4D',a*a-4*D,Q(9,10),Q(11,10),0,(4,8)),
 ('OOO P lower',P,Q(9,10),Q(997,1000),Q(1,10**9),(3,6)),
 ('OOO P upper',P,Q(1003,1000),Q(11,10),Q(1,10**9),(3,6))]: bern(name,expr,l,h,do,t,deg)
for j,sgn,target in [(0,1,Q(2,10**10)),(0,-1,Q(2,10**10)),(1,1,Q(3,10**11)),(1,-1,Q(7,10**12)),(2,1,Q(8,10**10))]:
 H=v*(sgn*r*D-n[j])**2+D*(D+a*a)*P
 bern('OOO H%d,%+d'%(j+1,sgn),H,Q(997,1000),Q(1003,1000),do,target,(9,20))
# Exact unshifted OOO endpoint margins and shift budgets.
bo=(1-v)*(v-Q(1,4))*(v-Q(1,2))/10
brackets=[Q(1,5),Q(3,10),Q(61,125),Q(64,125),Q(9,10),Q(11,10)]
expected=[Q(11,10000),Q(47,80000),Q(799,78125000),Q(9017,625000000),Q(193,80000),Q(1949,400000)]
for vi,margin in zip(brackets,expected):
    actual=abs(bo.subs(v,vi))-r*(vi/5+Q(3,40)*abs(1-vi))
    ck('OOO exact bracket margin '+str(vi),actual==margin and margin-32*do>Q(2,10**6))
ck('OOO first HB shift margin',-Q(15941,10**7)+64*do<-Q(3,2000))
ck('OOO second HB shift margin',Q(560777,7812500000)-64*do>Q(1,20000))

# EOO exact tensor Bernstein certificate
print('EOO table',flush=True)
a=d**4-3*d**3/2-6*d*d*v-d*d/4+9*d*v/2+d/2+v*v+v/4-1
b=-4*d**3+9*d*d/2+4*d*v+d/2-3*v/2-Q(1,2)
u=-17*d**3/8+2*d*d+51*d*v/8-15*d/8-2*v
dd=51*d*d/8-4*d-17*v/8+Q(15,8)
c=3*d**3/2-9*d*v/2+3*d/2; y=-9*d*d/2+3*v/2-Q(3,2)
C=v*dd*dd-c*c; P=C-a*a; n2=a*b*c-C*y
B=a*b+dd*u+c*y
Ab=a*a+u*u+c*c-v*(b*b+dd*dd+y*y)
for name,expr,l,h,t in [
 ('EOO C-1/20',C-Q(1,20),Q(99,100),Q(101,100),0),
 ('EOO 3/40-C',Q(3,40)-C,Q(99,100),Q(101,100),0),
 ('EOO 3/40-d2',Q(3,40)-dd*dd,Q(99,100),Q(101,100),0),
 ('EOO P lower',P,Q(99,100),Q(4997,5000),Q(3,2*10**9)),
 ('EOO P upper',P,Q(5003,5000),Q(101,100),Q(3,2*10**9))]: bern(name,expr,l,h,de,t)
for sig in [-1,1]:
 for tau in [-1,1]: bern('EOO face B %d %d'%(sig,tau),-sig*(B+dd*sig*r+c*tau*r),Q(4997,5000),Q(5003,5000),de,Q(17,10**5))
 H=C*(sig*r*C-n2)**2+dd*dd*C*(C+v*b*b)*P
 bern('EOO H2,%+d'%sig,H,Q(4997,5000),Q(5003,5000),de,Q(7,10**13),(12,26))
H3=v*r*r*C+(C+v*b*b)*P
bern('EOO H3',H3,Q(4997,5000),Q(5003,5000),de,Q(4,10**8),(7,14))
# EOO endpoint signs and first two interval auxiliary bounds
for vi,sgn in zip([Q(3,8),Q(2,5),Q(17,20),Q(7,8),Q(99,100),Q(101,100)],[1,-1,-1,1,1,-1]):
 for t1 in [-r,r]:
  for t2 in [-r,r]: bern('EOO bracket B %s %s %s'%(vi,t1,t2),(sgn*(B+dd*t1+c*t2)).subs(v,vi),0,1,de,Q(9,10000))
for lo,hi,a0,umax,vymax in [(Q(3,8),Q(2,5),-Q(1,25),1,Q(2,5)),(Q(17,20),Q(7,8),Q(3,20),2,Q(3,10))]:
 bern('EOO Ab %s'%lo, a0-Ab if a0<0 else Ab-a0,lo,hi,de)
 for sig in [-1,1]:
  bern('EOO u %s %d'%(lo,sig),umax-sig*u,lo,hi,de)
  bern('EOO vy %s %d'%(lo,sig),vymax-sig*v*y,lo,hi,de)
# analytic fixtures
alpha,beta=S.symbols('alpha beta',real=True)
N=S.expand((s+1)**4+alpha**2*s*s+beta**2)
ck('running norm',S.expand(norm([(s+1)**2,alpha*s,beta,0]).as_expr()-N)==0)
ck('one-component determinants',S.expand(4*(2*alpha**4+alpha**3-7*alpha**2+5*alpha-1)).subs(alpha,Q(9,20))==Q(451,20000))
ck('one-component unstable interior determinant',S.expand(4*(2*alpha**4+alpha**3-7*alpha**2+5*alpha-1)).subs(alpha,Q(3,4))==-Q(17,32))
B=-Q(3,4); C=-Q(13,50); dc=1-C*C; U=-B*C; V=B
Ls=dc**2-V**2-C*C*v; Fs=S.expand(Ls**2-(U*(dc-V)-C*v)**2-v*(dc+V-C*U)**2)
ck('Schur cubic L endpoint',Ls.subs(v,Q(3,4))==Q(1601061,6250000))
ck('Schur cubic F endpoints',Fs.subs(v,Q(1,4))==Q(1647104154759,156250000000000) and Fs.subs(v,Q(3,4))==Q(5454247783509,156250000000000),{'L':str(Ls),'F':str(Fs)})
rho=S.symbols('rho'); Fc=5-Q(11,2)*rho+Q(23,25)*rho*rho-Q(51,1000)*rho**3
ck('application endpoint gaps',Fc.subs(rho,1)==Q(369,1000) and Fc.subs(rho,Q(11,10))==-Q(4681,10**6))
ck('application isolator',Fc.subs(rho,Q(1098721768495,10**12))>0 and Fc.subs(rho,Q(1098721768496,10**12))<0)

# Exact root-count/isolation checks for both reserve polynomials.
ck('application unique reserve isolator',S.Poly(51*rho**3-920*rho**2+5500*rho-5000,rho).count_roots(Q(1098721768495,10**12),Q(1098721768496,10**12))==1)
ck('one-component reserve isolator',S.Poly(rho**3-30*rho**2+1100*rho-5000,rho).count_roots(Q(5143205348290,10**12),Q(5143205348291,10**12))==1)
ck('one-component u isolator',S.Poly(rho**3+8*rho-4,rho).count_roots(Q(4,9),Q(1,2))==1)
ck('one-component reserve substitution',S.expand((10*(1-rho))**3-30*(10*(1-rho))**2+1100*(10*(1-rho))-5000)==S.expand(-1000*(rho**3+8*rho-4)))
ck('one-component design gap',Q(9,10)*(Q(9,10)**2+8)-4==Q(3929,1000))
out['all_exact_checks_pass']=all(x['pass'] for x in out['exact_checks']) and all(x['pass'] for x in out['bernstein'])
if not out['all_exact_checks_pass']: raise AssertionError('An exact algebraic claim failed')
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--build',action='store_true',help='regenerate retained rational certificate')
ap.add_argument('--certificate',type=pathlib.Path,default=root/'exact_example_certificate.json')
args=ap.parse_args()
if args.build:
    args.certificate.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
else:
    saved=json.loads(args.certificate.read_text(encoding='utf-8'))
    if saved!=out: raise AssertionError('Recomputed exact data differ from retained certificate')
print('PASS exact algebraic checks:',len(out['exact_checks']))
print('PASS exact tensor Bernstein inequalities:',len(out['bernstein']))
