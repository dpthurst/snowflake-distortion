"""Approximate similarities of shortenings to the full arc and S configurations."""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ALPHA = np.log(3)/np.log(4)
FILES = {
    'S_plain': ROOT/'.ultramath/smwhdqmht/m1/s_m1600_plain.npy',
    'S_kick': ROOT/'.ultramath/smwhdqmht/m1/s_m1600_kick.npy',
    'arc': ROOT/'.ultramath/smwhdqmht/m1/arc_m1600.npy',
    'point_symmetric_200': ROOT/'.ultramath/smwhdqmht/m2/w_smoothS_t2_200.npy',
    'point_symmetric_400': ROOT/'.ultramath/smwhdqmht/m2/w_smoothS_t2_400.npy',
}

def gauge(v):
    v = np.asarray(v, dtype=complex)
    return (v-v[0])/(v[-1]-v[0])

def interp(v, n):
    v = gauge(v)
    t = np.linspace(0, 1, n)
    x = np.linspace(0, 1, len(v))
    return np.interp(t, x, v.real) + 1j*np.interp(t, x, v.imag)

def geometry(v):
    v = gauge(v); m = len(v)
    i, j = np.triu_indices(m, 1)
    r = np.abs(v[i]-v[j])/(j-i)**ALPHA
    dia = np.abs(v[:, None]-v[None, :]).max()
    return dict(D=float(r.max()/r.min()), rmin=float(r.min()), rmax=float(r.max()),
                diameter=float(dia), minsep=float(np.abs(v[:, None]-v[None, :]+np.eye(m)*1e9).min()),
                im_min=float(v.imag.min()), im_max=float(v.imag.max()),
                apex=float(v.imag.argmax()/(m-1)), trough=float(v.imag.argmin()/(m-1)))

def inverse_deviation(v):
    v = gauge(v); dia=np.abs(v[:, None]-v[None, :]).max()
    e = np.abs(v+v[::-1]-1)
    root = np.sqrt(np.mean(e**2))
    return dict(sup=float(e.max()/dia), rms=float(root/dia))

def mirror_deviation(v):
    v = gauge(v); dia=np.abs(v[:,None]-v[None,:]).max()
    e = np.abs(v-(1-np.conj(v[::-1])))
    return dict(sup=float(e.max()/dia), rms=float(np.sqrt(np.mean(e**2))/dia))

def proc(full, block, variant):
    if variant & 1:
        block = np.conj(block)
    if variant & 2:
        block = block[::-1]
    a, b = full-full.mean(), block-block.mean()
    coeff = np.vdot(a, b)/np.vdot(a, a)
    shift = block.mean()-coeff*full.mean()
    residual = block-(coeff*full+shift)
    dia = np.abs(np.subtract.outer(block, block)).max()
    return dict(coeff_real=float(coeff.real), coeff_imag=float(coeff.imag),
                max_err=float(np.abs(residual).max()/dia),
                rms_err=float(np.sqrt(np.mean(np.abs(residual)**2))/dia))

def block_stats(full_v, a, b, ncoarse=64):
    ix = np.rint(np.linspace(a, b, ncoarse)).astype(int)
    block = full_v[ix]
    fir = interp(full_v, ncoarse)
    out=[dict(kind=kind, **proc(fir, block, kind)) for kind in range(4)]
    unique_ix = ix[np.r_[True, np.diff(ix)>0]]
    i, j = np.triu_indices(len(unique_ix), 1)
    r = np.abs(full_v[unique_ix[i]]-full_v[unique_ix[j]])/(unique_ix[j]-unique_ix[i])**ALPHA
    return dict(interval=[float(ix[0]/(len(full_v)-1)), float(ix[-1]/(len(full_v)-1))],
                nodes=[int(ix[0]), int(ix[-1])], block_D=float(r.max()/r.min()), fits=out)

def choose(row):
    return min(row['fits'], key=lambda q:q['max_err'])['kind']

def sort_records(records):
    return sorted(records, key=lambda r: min(q['max_err'] for q in r['fits']))

def refine_matches(full_v, records, ncoarse=512):
    top = records[:12]
    for row in top:
        fresh = block_stats(full_v, row['nodes'][0], row['nodes'][1], ncoarse)
        row['refined_512'] = fresh['fits']
        row['block_D_refined'] = fresh['block_D']
    return top

def main():
    output={}
    long_lines=[]
    for name, path in FILES.items():
        v=gauge(np.load(path))
        g=geometry(v)
        print('\n===',name,'===\n',g)
        print('inversion (T2) deviation:',inverse_deviation(v),
              'mirror (T1) deviation:',mirror_deviation(v))
        all_blocks=[]
        Ks=[2,4,8,16,32]
        for K in Ks:
            for q in range(K):
                a=round(q*(len(v)-1)/K); b=round((q+1)*(len(v)-1)/K)
                row=block_stats(v,a,b)
                row['relative_length']=1/K
                row['kind']=choose(row)
                all_blocks.append(row)
        all_blocks=sort_records(all_blocks)
        by_length={str(k):[dict(interval=r['interval'],kind=r['kind'],block_D=r['block_D'],
                                max_err=r['fits'][r['kind']]['max_err'],rms_err=r['fits'][r['kind']]['rms_err'])
                           for r in all_blocks if r['relative_length']==1/k] for k in Ks}
        for k in Ks:
            z=by_length[str(k)]
            print('K',k,'best,median,max of max_err:',z[0]['max_err'],np.median([x['max_err'] for x in z]),z[-1]['max_err'])
            print('  best interval/type:',z[0]['interval'],z[0]['kind'],'similarity range sklearn[best three]:',
                  [(x['interval'],x['kind'],round(x['max_err'],5)) for x in z[:3]])
        top=refine_matches(v,all_blocks)
        for row in top[:6]:
            print(' refined 512',row['relative_length'],row['interval'],'kind',row['kind'],
                  'maxerr',row['refined_512'][row['kind']]['max_err'],
                  'rmserr',row['refined_512'][row['kind']]['rms_err'],'blockD',row['block_D_refined'])
        output[name]=dict(geometry=g,inversion=inverse_deviation(v),mirror=mirror_deviation(v),
                          by_length=by_length,best_512=top)
    (ROOT/'results/shortening-similarity.json').write_text(json.dumps(output,indent=2)+'\n')
if __name__=='__main__':main()
