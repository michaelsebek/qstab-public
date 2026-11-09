#!/usr/bin/env python3
"""Exact checks for the four-mode gyroscopic/circulatory PD benchmark.

Python standard library only. No simulation, eigenvalue sampling or measured
data are used. The proof of uniform positivity is in mechanical_details.tex;
this script verifies its polynomial identities and rational inequalities.
"""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import json
import sys

if not __debug__:
    raise RuntimeError('Run without -O: exact assertions are part of the checker.')

checks = []
def check(name, value):
    assert value, name
    checks.append(name)

def parity(p):
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1,len(p)))

def det(A):
    total = F(0)
    for p in permutations(range(len(A))):
        term = F(parity(p))
        for i,j in enumerate(p):
            term *= A[i][j]
        total += term
    return total

def matmul(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]

def left(a):
    r,i,j,k = a
    return [[r,-i,-j,-k],[i,r,-k,j],[j,k,r,-i],[k,-j,i,r]]

Li,Lj,Lk = [left(v) for v in [(0,1,0,0),(0,0,1,0),(0,0,0,1)]]
I = [[int(i==j) for j in range(4)] for i in range(4)]
for name,L in zip(('i','j','k'),(Li,Lj,Lk)):
    check('skew_'+name, all(L[i][j] == -L[j][i] for i in range(4) for j in range(4)))
    check('square_'+name, matmul(L,L) == [[-z for z in row] for row in I])
check('quaternion_product_ij', matmul(Li,Lj) == Lk)
check('quaternion_product_ji', matmul(Lj,Li) == [[-z for z in row] for row in Lk])

# Sparse exact polynomials in (d,b,u,v); indices are exponent tuples.
ZERO=(0,0,0,0)
def P(x): return {} if x == 0 else {ZERO:F(x)}
def add(a,b):
    out = dict(a)
    for e,c in b.items(): out[e] = out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}
def neg(a): return {e:-c for e,c in a.items()}
def sub(a,b): return add(a,neg(b))
def mul(a,b):
    out={}
    for e,c in a.items():
        for f,z in b.items():
            t=tuple(x+y for x,y in zip(e,f))
            out[t]=out.get(t,F(0))+c*z
    return {e:c for e,c in out.items() if c}
def scale(a,x): return mul(a,P(x))
def power(a,n):
    out=P(1)
    for _ in range(n): out=mul(out,a)
    return out
def variable(j):
    e=[0]*4; e[j]=1
    return {tuple(e):F(1)}
def peval(a,vals):
    return sum(c*product(vals[k]**e[k] for k in range(4)) for e,c in a.items())
def product(seq):
    ans=F(1)
    for x in seq: ans*=x
    return ans
def pdet(A):
    total={}
    for p in permutations(range(len(A))):
        term=P(parity(p))
        for i,j in enumerate(p): term=mul(term,A[i][j])
        total=add(total,term)
    return total

d,b,u,v=[variable(j) for j in range(4)]
r2=add(add(power(b,2),power(u,2)),power(v,2))
coeff=[P(1),scale(d,2),add(add(power(d,2),P(4)),r2),scale(add(d,b),4),P(8)]
H=[[coeff[1],coeff[3],{},{}],
   [P(1),coeff[2],coeff[4],{}],
   [{},coeff[1],coeff[3],{}],
   [{},P(1),coeff[2],coeff[4]]]
minors=[pdet([row[:j] for row in H[:j]]) for j in range(1,5)]
A=add(add(scale(power(d,2),2),scale(mul(d,b),2)),P(-4))
expected_D2=scale(sub(add(add(power(d,3),scale(d,2)),mul(d,r2)),scale(b,2)),2)
expected_D3=scale(add(mul(add(power(d,2),power(b,2)),A),
                         scale(mul(mul(d,add(d,b)),add(power(u,2),power(v,2))),2)),4)
check('Delta1_identity', minors[0] == scale(d,2))
check('Delta2_identity', minors[1] == expected_D2)
check('Delta3_identity', minors[2] == expected_D3)
check('Delta4_identity', minors[3] == scale(expected_D3,8))

verts=[[F(x,10) for x in row] for row in
       [[3,5,0,0],[2,0,5,0],[2,0,0,5],[2,-5,-5,-5],[1,0,0,0]]]
diff=[[verts[j][i]-verts[4][i] for j in range(4)] for i in range(4)]
rank_minor=det(diff)
check('affine_dimension_4', abs(rank_minor)==F(1,16))
centroid=[sum(row[i] for row in verts)/5 for i in range(4)]
check('centroid', centroid==[F(1,5),F(0),F(0),F(0)])

vertex_minors=[]
for n,row in enumerate(verts,1):
    c,gi,gj,gk=row
    check('c_bounds_vertex_'+str(n),F(1,10)<=c<=F(3,10))
    check('g1_bounds_vertex_'+str(n),F(-1,2)<=gi<=F(1,2))
    check('correlation_vertex_'+str(n),gi+5*(c-F(1,10))>=0)
    vals=[c+F(3,2),gi,gj,gk]
    dm=[peval(p,vals) for p in minors]
    # Independent numerical-coefficient determinant for every vertex.
    dd,bb,uu,vv=vals
    co=[F(1),2*dd,dd*dd+4+bb*bb+uu*uu+vv*vv,4*(dd+bb),F(8)]
    Hn=[[co[1],co[3],0,0],[1,co[2],co[4],0],[0,co[1],co[3],0],[0,1,co[2],co[4]]]
    check('vertex_determinants_'+str(n), dm==[det([r[:j] for r in Hn[:j]]) for j in range(1,5)])
    check('vertex_Hurwitz_'+str(n),all(x>0 for x in dm))
    vertex_minors.append([str(x) for x in dm])

# Exact interval proof for A = 2d^2+2dg1-4.
# For c in [1/10,1/5], A >= 2-8c-8c^2; its derivative is negative.
# For c in [1/5,3/10], A >= 2c^2+5c-1; its derivative is positive.
cL,cM,cR=F(1,10),F(1,5),F(3,10)
left=lambda c:2-8*c-8*c*c
right=lambda c:2*c*c+5*c-1
check('left_decreasing',-8-16*cL<0 and -8-16*cR<0)
check('right_increasing',4*cL+5>0 and 4*cR+5>0)
check('piecewise_A_lower_bound',left(cM)==right(cM)==F(2,25))
dmin=F(8,5)
bmin=F(-1,2)
check('positive_second_Delta3_term',dmin>0 and dmin+bmin>0)
lower=[2*dmin,2*(dmin**3+2*dmin-1),4*dmin*dmin*F(2,25)]
lower.append(8*lower[2])
check('uniform_Hurwitz_lower_bounds',lower==[F(16,5),F(1574,125),F(512,625),F(4096,625)])

# The uncontrolled centroid has k=1, h=2, d=c=1/5, g=0.
oc=[F(1),F(2,5),F(51,25),F(2,5),F(5)]
Ho=[[oc[1],oc[3],0],[1,oc[2],oc[4]],[0,oc[1],oc[3]]]
open_D3=det(Ho)
check('uncontrolled_centroid_unstable',open_D3==F(-396,625))

# This box point is inadmissible in the correlated simplex (minimum c forces g=0).
box=[F(8,5),F(-1,2),F(0),F(0)]
box_D3=peval(minors[2],box)
check('independent_box_unstable',box_D3==F(-3372,625))
check('box_witness_not_admissible',F(-1,2)+5*(F(1,10)-F(1,10))<0)

report={
 'status':'ALL EXACT MECHANICAL CHECKS OK',
 'check_count':len(checks),
 'checks':checks,
 'normalized_parameters':{'mass':1,'potential_stiffness':1,'circulatory_h':['0','2','0','0'],
   'PD_position_gain':'1','PD_velocity_gain':'3/2'},
 'uncertainty_vertices_c_g1_g2_g3':[[str(z) for z in row] for row in verts],
 'affine_rank_minor':str(rank_minor),
 'affine_dimension':4,
 'uniform_Hurwitz_lower_bounds':[str(z) for z in lower],
 'vertex_Hurwitz_minors':vertex_minors,
 'uncontrolled_centroid_Delta3':str(open_D3),
 'independent_box_witness':{'c':'1/10','g':['-1/2','0','0'],'closed_loop_Delta3':str(box_D3)},
 'data_origin':'Explicitly constructed idealized mechanical benchmark; no measured or MATLAB data.',
 'verification':'Exact sparse polynomial identities and rational determinant/inequality checks; no sampled stability proof.'
}
target=Path(__file__).with_name('mechanical_checks.json')
target.write_text(json.dumps(report,indent=2)+'\n')
print(report['status'])
print('checks =',len(checks))
print('affine determinant =',rank_minor)
print('uniform Hurwitz lower bounds =',report['uniform_Hurwitz_lower_bounds'])
print('uncontrolled centroid Delta3 =',open_D3)
print('independent-box Delta3 =',box_D3)
