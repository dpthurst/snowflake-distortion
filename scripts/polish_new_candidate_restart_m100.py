"""Continue the one incomplete SLP polish and refresh the saved restart records."""
import json
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'.ultramath/smwgx79ya/m8'))
from core2 import Problem, norm_shape, supdist
from slp2 import slp, stress, grad_rows
from restart_new_candidate_m100 import stats

path=ROOT/'results/new_candidate_uniform_m100_restarts.npz'
arrays=dict(np.load(path))
records=json.loads((ROOT/'results/new_candidate_uniform_m100_restarts.json').read_text())
P=Problem(100,gauge='v2')
label='perturb_1e-5_seed_41'
r=next(x for x in records if x['label']==label)
v=arrays[label+'_polished']
for batch in range(4):
    v=(v-v[0])/(v[1]-v[0])
    v,_=slp(P,v,rho0=0.001,iters=500)
    v=norm_shape(v)
    st=stress(P,v,tol=1e-8)
    print('extra polish batch',batch+1,'D',stats(v)['D'],'stress',st['status'],
          'tight',len(st['pairs_max']),len(st['pairs_min']),flush=True)
    arrays[label+'_polished']=v
    if st['status']==0:break
A,L=grad_rows(P,v)
residual=None
if st['status']==0:
    residual=st['lam']@A[st['imax']].toarray()-st['mu']@A[st['imin']].toarray()
    arrays[label+'_multipliers_max']=st['lam']
    arrays[label+'_multipliers_min']=st['mu']
else:
    arrays.pop(label+'_multipliers_max',None)
    arrays.pop(label+'_multipliers_min',None)
r.update(polished=stats(v),tight_max=len(st['pairs_max']),tight_min=len(st['pairs_min']),
        stress_status=int(st['status']),
        multiplier_residual_inf=None if residual is None else float(np.max(np.abs(residual))),
        pairs_max=[list(map(int,p)) for p in st['pairs_max']],
        pairs_min=[list(map(int,p)) for p in st['pairs_min']])
cloud=np.load(ROOT/'.scratch/final_m100.npz',allow_pickle=True)
old=cloud['W'][np.argmin(cloud['Ds'])]
for record in records:
    label=record['label']
    record['old_best_supdist']=float(supdist(arrays[label+'_polished'],old))
    record['polished']=stats(arrays[label+'_polished'])
    print(label,'final',record['polished']['D'],'distance to old best',
          record['old_best_supdist'],flush=True)
np.savez(path,**arrays)
(ROOT/'results/new_candidate_uniform_m100_restarts.json').write_text(
    json.dumps(records,indent=2,allow_nan=False)+'\n')
