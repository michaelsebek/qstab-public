#!/usr/bin/env python3
"""Build an exact local certificate for a five-dimensional cubic simplex.

Python standard library only. The initial polarization hull contains zero at
frequency zero. Longest-edge parameter bisection and exact frequency Bernstein
separation certify the complete family. All proof arithmetic is rational.
"""
from __future__ import annotations
from fractions import Fraction as Q
from itertools import combinations
from math import comb, isqrt, gcd, lcm
from pathlib import Path
import argparse, json, time

ZERO=(Q(0),)*4
ONE=(Q(1),Q(0),Q(0),Q(0))

def qadd(a,b): return tuple(x+y for x,y in zip(a,b))
def qmul(a,b):
 w,x,y,z=a;v,r,s,t=b
 return (w*v-x*r-y*s-z*t,w*r+x*v+y*t-z*s,w*s-x*t+y*v+z*r,w*t+x*s-y*r+z*v)
def mul(a,b):
 c=[ZERO]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b): c[i+j]=qadd(c[i+j],qmul(x,y))
 return c

def add(a,b):
 return [qadd(a[i] if i<len(a) else ZERO,b[i] if i<len(b) else ZERO) for i in range(max(len(a),len(b)))]
def scale(a,t): return [tuple(t*x for x in c) for c in a]
def midpoint(a,b): return [tuple((x+y)/2 for x,y in zip(c,d)) for c,d in zip(a,b)]
def polarization(a,b):
 return tuple(sum((sum(a[r][nu]*b[k-r][nu] for nu in range(4)) for r in range(len(a)) if 0<=k-r<len(b)),Q(0)) for k in range(len(a)+len(b)-1))
def parts(c):
 return tuple(x*(-1)**(k//2) if k%2==0 else Q(0) for k,x in enumerate(c)),tuple(x*(-1)**((k-1)//2) if k%2 else Q(0) for k,x in enumerate(c))
def bernstein(c,a,b):
 n=len(c)-1;w=b-a
 t=[w**j*sum((c[k]*comb(k,j)*a**(k-j) for k in range(j,n+1)),Q(0)) for j in range(n+1)]
 return tuple(sum((t[j]*Q(comb(i,j),comb(n,j)) for j in range(i+1)),Q(0)) for i in range(n+1))
def value(c,w):
 ans=Q(0)
 for x in reversed(c):ans=ans*w+x
 return ans

def det(a):
 a=[list(r) for r in a];n=len(a);d=Q(1)
 for i in range(n):
  p=next((j for j in range(i,n) if a[j][i]),None)
  if p is None:return Q(0)
  if p!=i:a[i],a[p]=a[p],a[i];d=-d
  x=a[i][i];d*=x
  for j in range(i+1,n):
   f=a[j][i]/x
   for k in range(i+1,n):a[j][k]-=f*a[i][k]
   a[j][i]=0
 return d

def hurwitz(c):
 n=len(c)-1
 h=[[c[n+i-2*j] if 0<=n+i-2*j<=n else Q(0) for j in range(1,n+1)] for i in range(1,n+1)]
 return [det([r[:k] for r in h[:k]]) for k in range(1,n+1)]
def ceil_sqrt(x):
 a=isqrt(x.numerator//x.denominator)
 return a+int(a*a<x)
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v) for v in x]
 return x

def family():
 A=[(Q(1),Q(1),Q(0),Q(0)),ONE]
 B=[(Q(3),Q(0),Q(0),Q(1)),ONE]
 signs=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
 cs=[(Q(2),*(Q(3*s) for s in v)) for v in signs]
 ps=[mul(mul(A,[c,ONE]),B) for c in cs]
 ps.append(add(ps[0],scale(mul(A,B),Q(1,10))))
 ps.append(add(ps[0],[ZERO,ZERO,tuple(Q(1,10)*v for v in ONE)]))
 return A,B,cs,ps


def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
def separator(points):
 """Rational strict separator of a planar point hull, or None."""
 p=sorted(set(points))
 candidates=[(1,0),(-1,0),(0,1),(0,-1)]
 if len(p)==1:candidates.append(p[0])
 elif len(p)==2:
  a,b=p;v=(b[0]-a[0],b[1]-a[1]);vv=v[0]*v[0]+v[1]*v[1]
  t=max(Q(0),min(Q(1),-(a[0]*v[0]+a[1]*v[1])/vv))
  candidates.append((a[0]+t*v[0],a[1]+t*v[1]))
 else:
  lower=[];upper=[]
  for q in p:
   while len(lower)>1 and cross(lower[-2],lower[-1],q)<=0:lower.pop()
   lower.append(q)
  for q in reversed(p):
   while len(upper)>1 and cross(upper[-2],upper[-1],q)<=0:upper.pop()
   upper.append(q)
  hull=lower[:-1]+upper[:-1]
  for a,b in zip(hull,hull[1:]+hull[:1]):candidates.append((a[1]-b[1],b[0]-a[0]))
 viable=[]
 for x,y in candidates:
  size=abs(x)+abs(y)
  if not size:continue
  gap=min(x*a+y*b for a,b in p)
  if gap>0:viable.append((gap/size,(Q(x),Q(y))))
 if not viable:return None
 _,h=max(viable)
 # A small integer normal makes the serialized certificate compact.
 size=abs(h[0])+abs(h[1])
 for grid in [1,2,4,8,16,32,64,128,256,512,1024,4096,65536]:
  v=(round(grid*h[0]/size),round(grid*h[1]/size))
  if min(v[0]*a+v[1]*b for a,b in p)>0:return v
 d=lcm(h[0].denominator,h[1].denominator);x=int(h[0]*d);y=int(h[1]*d);g=gcd(abs(x),abs(y))
 return x//g,y//g


def main():
 if not __debug__:raise RuntimeError('Run without -O: exact assertions are required.')
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--out',default=str(Path(__file__).resolve().parent))
 parser.add_argument('--max-parameter-depth',type=int,default=16)
 parser.add_argument('--max-frequency-depth',type=int,default=18)
 args=parser.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
 A,B,cs,ps=family()
 minors=[hurwitz(polarization(p,p)) for p in ps]
 assert all(x>0 for row in minors for x in row)
 bounds2=[max(sum(x*x for x in p[k]) for p in ps) for k in range(3)]
 bounds=[ceil_sqrt(x) for x in bounds2];R=Q(1+max(bounds))
 flat=[[v for q in p[:3] for v in q] for p in ps]
 rank_minor=None
 for rows in combinations(range(12),5):
  d=det([[flat[j][r]-flat[5][r] for j in range(5)] for r in rows])
  if d:rank_minor=(rows,d);break
 assert rank_minor
 b11=polarization(ps[0],ps[0])[0];b12=polarization(ps[0],ps[1])[0]
 assert b11==620 and b12==-100 and Q(5,36)*b11+Q(31,36)*b12==0
 stats={'parameter_leaves':0,'parameter_splits':0,'frequency_cells':0,'maximum_parameter_depth':0,'maximum_frequency_depth':0,'bernstein_identities':0}
 global_gap=None
 control_cache={}
 def polarkey(a,b):
  aa=tuple(a);bb=tuple(b)
  return (aa,bb) if aa<=bb else (bb,aa)
 def controls(qs):
  ans=[]
  for i in range(4):
   for j in range(i,4):
    key=polarkey(qs[i],qs[j])
    if key not in control_cache:control_cache[key]=parts(polarization(qs[i],qs[j]))
    ans.append(control_cache[key])
  return ans
 def frequency_certificate(ctrl):
  pending=[(Q(0),R,0)];cells=[];max_depth=0
  while pending:
   a,b,depth=pending.pop();middle=(a+b)/2
   pts=[(value(re,middle),value(im,middle)) for re,im in ctrl]
   h=separator(pts)
   if h is None:return None
   scalar=[tuple(h[0]*x+h[1]*y for x,y in zip(re,im)) for re,im in ctrl]
   beta=[bernstein(c,a,b) for c in scalar]
   lower=min(min(row) for row in beta)
   if lower>0:
    gap=lower/(abs(h[0])+abs(h[1]))
    cells.append({'left':a,'right':b,'h':h,'bernstein':beta,'modulus_lower_bound':gap,'frequency_depth':depth})
    max_depth=max(max_depth,depth)
   else:
    if depth>=args.max_frequency_depth:return None
    pending.extend([(middle,b,depth+1),(a,middle,depth+1)])
  cells.sort(key=lambda c:c['left'])
  return cells,max_depth
 def recurse(qs,depth=0):
  nonlocal global_gap
  ctrl=controls(qs)
  # At zero all controls are real. Avoid futile frequency refinement if the
  # coarse cell has a genuine zero-containing control interval already there.
  if min(re[0] for re,im in ctrl)<=0:
   result=None
  else:result=frequency_certificate(ctrl)
  if result is not None:
   cells,fd=result;stats['parameter_leaves']+=1;stats['frequency_cells']+=len(cells);stats['bernstein_identities']+=10*len(cells)
   stats['maximum_parameter_depth']=max(stats['maximum_parameter_depth'],depth);stats['maximum_frequency_depth']=max(stats['maximum_frequency_depth'],fd)
   gap=min(c['modulus_lower_bound'] for c in cells);global_gap=gap if global_gap is None else min(global_gap,gap)
   return {'cells':cells}
  if depth>=args.max_parameter_depth:raise RuntimeError('Unresolved parameter simplex; no robustness claim.')
  # Longest edge in the complete real coefficient space, deterministic ties.
  i,j=max(combinations(range(4),2),key=lambda ij:sum((x-y)**2 for c,d in zip(qs[ij[0]],qs[ij[1]]) for x,y in zip(c,d)))
  m=midpoint(qs[i],qs[j]);left=list(qs);right=list(qs);left[i]=m;right[j]=m
  stats['parameter_splits']+=1
  return {'split_edge':[i,j],'children':[recurse(left,depth+1),recurse(right,depth+1)]}
 faces=[]
 for vertices in combinations(range(6),4):
  before=stats.copy();tree=recurse([ps[j] for j in vertices]);faces.append({'vertices':[j+1 for j in vertices],'tree':tree})
  print('face',tuple(j+1 for j in vertices),'leaves',stats['parameter_leaves']-before['parameter_leaves'],'cells',stats['frequency_cells']-before['frequency_cells'],flush=True)
 data={'example':'five-dimensional cubic simplex with a strict global-hull false alarm',
       'generators':ps,'outer_left':A,'outer_right':B,'middle_constants':cs,
       'perturbation':Q(1,10),'affine_minor_rows':rank_minor[0],'affine_minor':rank_minor[1],
       'vertex_hurwitz_minors':minors,'coefficient_norm_squared_maxima':bounds2,'coefficient_norm_upper_bounds':bounds,'frequency_cutoff':R,
       'global_false_alarm':{'frequency':Q(0),'controls':[[1,1],[1,2]],'values':[b11,b12],'weights':[Q(5,36),Q(31,36)]},
       'faces':faces,'statistics':stats,'uniform_modulus_lower_bound_on_tested_faces':global_gap,
       'status':'CERTIFIED ROBUST HURWITZ WITH PARAMETER SUBDIVISION'}
 (out/'local_certificate.json').write_text(json.dumps(enc(data),indent=2)+'\n')
 summary={k:v for k,v in data.items() if k in ['example','statistics','uniform_modulus_lower_bound_on_tested_faces','frequency_cutoff','affine_minor','status']}
 (out/'local_certificate_summary.json').write_text(json.dumps(enc(summary),indent=2)+'\n')
 print(json.dumps(enc(summary),indent=2));print('Informational elapsed seconds:',time.perf_counter()-start)
if __name__=='__main__':main()
