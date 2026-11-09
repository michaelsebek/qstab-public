#!/usr/bin/env python3
"""Independent exact replay of local_certificate.json.

This program does not import the generator. It reconstructs the specified
quaternion polynomials using left multiplication matrices, verifies all
parameter bisections recursively, and expands every recorded frequency
Bernstein representation back to its exact power polynomial. Standard library
only; no frequency or parameter sampling is part of the proof.
"""
from fractions import Fraction as Q
from itertools import combinations,permutations
from math import comb
from pathlib import Path
import json,sys,time


def qproduct(a,b):
 w,x,y,z=a
 L=((w,-x,-y,-z),(x,w,-z,y),(y,z,w,-x),(z,-y,x,w))
 return tuple(sum(L[i][j]*b[j] for j in range(4)) for i in range(4))

def pproduct(a,b):
 z=[(Q(0),)*4 for _ in range(len(a)+len(b)-1)]
 for k in range(len(z)):
  z[k]=tuple(sum((qproduct(a[r],b[k-r])[j] for r in range(len(a)) if 0<=k-r<len(b)),Q(0)) for j in range(4))
 return z

def det(a):
 n=len(a);ans=Q(0)
 for p in permutations(range(n)):
  term=Q((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))
  for i,j in enumerate(p):
   term*=a[i][j]
   if not term:break
  ans+=term
 return ans

def conv(a,b):
 z=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):z[i+j]+=x*y
 return z

def polarpower(p,q):
 z=[Q(0)]*7
 for j in range(4):
  v=conv([c[j] for c in p],[c[j] for c in q]);z=[x+y for x,y in zip(z,v)]
 return z

def polynomial_power(p,n):
 z=[Q(1)]
 for _ in range(n):z=conv(z,p)
 return z

def basis_matrix(a,b):
 return [[x*Q(comb(6,i),(b-a)**6) for x in conv(polynomial_power([-a,Q(1)],i),polynomial_power([b,Q(-1)],6-i))] for i in range(7)]

def reconstruct_family():
 e=(Q(1),Q(0),Q(0),Q(0));z=(Q(0),)*4
 A=[(Q(1),Q(1),Q(0),Q(0)),e];B=[(Q(3),Q(0),Q(0),Q(1)),e]
 cs=[(Q(2),Q(3),Q(3),Q(3)),(Q(2),Q(3),Q(-3),Q(-3)),(Q(2),Q(-3),Q(3),Q(-3)),(Q(2),Q(-3),Q(-3),Q(3))]
 ps=[pproduct(A,pproduct([c,e],B)) for c in cs]
 ab=pproduct(A,B)
 ps.append([tuple(ps[0][k][j]+(ab[k][j]/10 if k<3 else 0) for j in range(4)) for k in range(4)])
 ps.append([tuple(ps[0][k][j]+(Q(1,10) if k==2 and j==0 else 0) for j in range(4)) for k in range(4)])
 return ps

def main():
 if not __debug__:raise RuntimeError('Run without -O: exact assertions are required.')
 path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('local_certificate.json')
 raw=json.loads(path.read_text());start=time.perf_counter()
 p=[[tuple(Q(x) for x in c) for c in q] for q in raw['generators']]
 assert len(p)==6 and all(len(q)==4 and all(len(c)==4 for c in q) for q in p)
 assert p==reconstruct_family()
 assert all(q[3]==(Q(1),Q(0),Q(0),Q(0)) for q in p)
 rows=raw['affine_minor_rows'];assert len(rows)==5 and len(set(rows))==5 and all(type(r)is int and 0<=r<12 for r in rows)
 flat=[[x for c in q[:3] for x in c] for q in p]
 affine_minor=det([[flat[j][r]-flat[5][r] for j in range(5)] for r in rows])
 assert affine_minor==Q(raw['affine_minor'])!=0
 global_controls={(i,j):polarpower(p[i],p[j]) for i in range(6) for j in range(i,6)}
 all_minors=[]
 for i in range(6):
  c=global_controls[(i,i)]
  h=[[c[6+r-2*s] if 0<=6+r-2*s<=6 else Q(0) for s in range(1,7)] for r in range(1,7)]
  minors=[det([r[:k] for r in h[:k]]) for k in range(1,7)]
  assert all(x>0 for x in minors);all_minors.append(minors)
 assert all_minors==[[Q(x) for x in row] for row in raw['vertex_hurwitz_minors']]
 norms2=[max(sum(x*x for x in q[k]) for q in p) for k in range(3)]
 assert norms2==[Q(x) for x in raw['coefficient_norm_squared_maxima']]
 bounds=raw['coefficient_norm_upper_bounds']
 assert len(bounds)==3 and all(type(x)is int and x>=0 for x in bounds)
 assert all(x*x>=v for x,v in zip(bounds,norms2))
 R=Q(raw['frequency_cutoff']);assert R==1+max(bounds)==30
 witness=raw['global_false_alarm'];assert Q(witness['frequency'])==0
 assert witness['controls']==[[1,1],[1,2]]
 values=[global_controls[(0,0)][0],global_controls[(0,1)][0]]
 weights=[Q(x) for x in witness['weights']]
 assert values==[Q(x) for x in witness['values']]==[Q(620),Q(-100)]
 assert weights==[Q(5,36),Q(31,36)] and all(w>0 for w in weights)
 assert sum(weights)==1 and sum(w*z for w,z in zip(weights,values))==0
 left_inverse=(Q(1,2),Q(-1,2),Q(0),Q(0))
 right_inverse=(Q(3,10),Q(0),Q(0),Q(-1,10))
 recovered=[qproduct(left_inverse,qproduct(q[0],right_inverse)) for q in p]
 assert min(q[0] for q in recovered)==2
 assert tuple(sum(q[j] for q in recovered[:4])/4 for j in range(4))==(Q(2),Q(0),Q(0),Q(0))
 actual_zero_bound=Q(20)*min(q[0] for q in recovered)**2
 assert actual_zero_bound==80
 stats={'parameter_leaves':0,'parameter_splits':0,'frequency_cells':0,'maximum_parameter_depth':0,'maximum_frequency_depth':0,'bernstein_identities':0}
 min_gap=None;min_zero_control=None;matrices={}
 def visit(node,qs,depth):
  nonlocal min_gap,min_zero_control
  if 'split_edge' in node:
   assert set(node)=={'split_edge','children'}
   edge=node['split_edge'];assert len(edge)==2 and all(type(x)is int and 0<=x<4 for x in edge) and edge[0]!=edge[1]
   assert len(node['children'])==2
   i,j=edge;m=[tuple((x+y)/2 for x,y in zip(a,b)) for a,b in zip(qs[i],qs[j])]
   left=list(qs);right=list(qs);left[i]=m;right[j]=m
   stats['parameter_splits']+=1
   # Replacing each endpoint in turn by the edge midpoint covers exactly the
   # original tetrahedron, with only the dividing face shared.
   visit(node['children'][0],left,depth+1);visit(node['children'][1],right,depth+1)
   return
  assert set(node)=={'cells'} and node['cells']
  stats['parameter_leaves']+=1;stats['maximum_parameter_depth']=max(stats['maximum_parameter_depth'],depth)
  ctr=[polarpower(qs[i],qs[j]) for i in range(4) for j in range(i,4)]
  assert len(ctr)==10
  local_zero=min(c[0] for c in ctr)
  min_zero_control=local_zero if min_zero_control is None else min(min_zero_control,local_zero)
  endpoint=Q(0)
  for cell in node['cells']:
   assert set(cell)=={'left','right','h','bernstein','modulus_lower_bound','frequency_depth'}
   a,b=Q(cell['left']),Q(cell['right']);assert a==endpoint and 0<=a<b<=R;endpoint=b
   fd=cell['frequency_depth'];assert type(fd)is int and fd>=0 and b-a==R/2**fd and (a/R*2**fd).denominator==1
   h=cell['h'];assert len(h)==2 and all(type(x)is int for x in h) and abs(h[0])+abs(h[1])>0
   beta=[[Q(x) for x in row] for row in cell['bernstein']]
   assert len(beta)==10 and all(len(row)==7 for row in beta)
   if (a,b) not in matrices:matrices[(a,b)]=basis_matrix(a,b)
   basis=matrices[(a,b)]
   for c,beta_row in zip(ctr,beta):
    assert min(beta_row)>0
    expanded=[sum((beta_row[i]*basis[i][k] for i in range(7)),Q(0)) for k in range(7)]
    expected=[c[k]*((-1)**(k//2))*h[0] if k%2==0 else c[k]*((-1)**((k-1)//2))*h[1] for k in range(7)]
    assert expanded==expected
    stats['bernstein_identities']+=1
   gap=min(x for row in beta for x in row)/(abs(h[0])+abs(h[1]))
   assert gap==Q(cell['modulus_lower_bound'])>0
   min_gap=gap if min_gap is None else min(min_gap,gap)
   stats['frequency_cells']+=1;stats['maximum_frequency_depth']=max(stats['maximum_frequency_depth'],fd)
  assert endpoint==R
 faces=set()
 for face in raw['faces']:
  assert set(face)=={'vertices','tree'}
  v=tuple(face['vertices']);assert v==tuple(sorted(v)) and len(v)==4 and len(set(v))==4 and all(1<=i<=6 for i in v) and v not in faces
  faces.add(v);visit(face['tree'],[p[j-1] for j in v],0)
 assert faces==set(combinations(range(1,7),4))
 assert stats==raw['statistics']
 assert stats['parameter_leaves']==stats['parameter_splits']+len(faces)
 assert min_gap==Q(raw['uniform_modulus_lower_bound_on_tested_faces'])>0
 result={'status':'INDEPENDENT EXACT LOCAL REPLAY OK','tetrahedral_subfamilies':len(faces),'positive_vertex_hurwitz_minors':36,'affine_minor':str(affine_minor),'frequency_cutoff':str(R),'global_hull_zero_witness_verified':True,'actual_zero_frequency_norm_lower_bound':str(actual_zero_bound),**stats,'uniform_modulus_lower_bound_on_tested_faces':str(min_gap),'distinct_frequency_intervals':len(matrices),'minimum_local_zero_frequency_control':str(min_zero_control),'global_zero_frequency_control_range':[str(min(c[0] for c in global_controls.values())),str(max(c[0] for c in global_controls.values()))]}
 path.with_name('local_certificate_replay.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2));print('Informational replay seconds:',time.perf_counter()-start)
if __name__=='__main__':main()
