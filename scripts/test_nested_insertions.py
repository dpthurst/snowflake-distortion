"""Test nested, parameter-compatible snowflake configurations by free insertion.
Unlike canonical two-substep refinement, inserted vertices are freely optimised;
old vertices stay fixed. Uses m=21 initially from the stored m=100 best shape.
"""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'.ultramath/smwgx79ya/m8'))
from core2 import Problem
A=np.log(3)/np.log(4)
def stats(v):
 i,j=np.triu_indices(len(v),1)
 r=abs(v[i]-v[j])/(j-i)**A
 return {'D':float(r.max()/r.min()),'rmin':float(r.min()),'rmax':float(r.max()),
         'min_step':float(abs(np.diff(v)).min()),'min_pair':float(abs(v[i]-v[j]).min())}
def insertion(v):
 m=len(v);w=np.empty(2*m-1,complex);w[::2]=v;w[1::2]=(v[:-1]+v[1:])/2
 return w

def new_level(previous,name):
 v=insertion(previous);n=len(v)
 # free_idx must contain only new (odd-index) vertices. Gauge 'vm' fixes 0 and 1.
 P=Problem(n,gauge='vm');old=P.free_idx.copy();P.free_idx=np.arange(1,n-1,2)
 P.nfree=2*len(P.free_idx);P.pos=-np.ones(n,int);P.pos[P.free_idx]=np.arange(len(P.free_idx))
 fixed=np.arange(0,n-1,2)
 assert np.all(P.pos[P.free_idx]>=0) and np.all(P.pos[fixed]<0)
 fixed_values=v.copy()
 def unpack_fixed(x):
  out=fixed_values.copy();n=len(P.free_idx)
  out[P.free_idx]=x[:n]+1j*x[n:]
  return out
 P.unpack=unpack_fixed
 x=P.pack(v)
 restored=P.unpack(x)
 assert np.max(abs(restored-v)) == 0, (name,P.free_idx[:5],P.free_idx[-5:],
                                         restored[:5],v[:5],restored[-5:],v[-5:])
 best=stats(restored)['D'];best_x=x.copy()
 print(name,'initial',stats(restored),flush=True)
 for beta in [0.1,0.3,1,3,10,30,100,300,1000,3000,10000,30000,100000]:
  for _ in range(2):
   r=minimize(P.fg,x,args=(beta,),jac=True,method='L-BFGS-B',
       options={'maxiter':350,'gtol':1e-14,'ftol':1e-16,'maxcor':30});x=r.x
   out=P.unpack(x);D=stats(out)['D']
   if np.isfinite(D) and D<best:best,best_x=D,x.copy()
  print(name,'beta',beta,'current',D,'best',best,flush=True)
  if not np.isfinite(D):x=best_x.copy();break
 v=P.unpack(best_x)
 # Exact preservation check on old vertices; objective determined independently.
 print(name,'FINAL',stats(v),'old vertex sup error',np.max(abs(v[::2]-previous)),flush=True)
 assert np.max(abs(v[::2]-previous))==0
 return v

def main():
 z=np.load(ROOT/'.scratch/final_m100.npz');w=z['W'][z['Ds'].argmin()]
 ids=np.rint(np.linspace(0,99,21)).astype(int);v=w[ids]
 v=(v-v[0])/(v[-1]-v[0])
 print('coarse start',stats(v),flush=True)
 arrays={'m21':v};records=[{'m':21,**stats(v)}]
 for m in (41,81,161):
  v=new_level(v,'m'+str(m));arrays['m'+str(m)]=v;records.append({'m':m,**stats(v)})
  np.savez(ROOT/'results/nested_insertion_tower.npz',**arrays)
  (ROOT/'results/nested_insertion_tower.json').write_text(json.dumps(records,indent=2)+'\n')
 print('D sequence',[(r['m'],r['D']) for r in records],flush=True)
if __name__=='__main__':main()
