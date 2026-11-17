import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
al=np.log(3)/np.log(4)
res=json.load(open('m40_clusters.json'))
def norm_end(V):
    W=(V-V[0])/(V[-1]-V[0])
    if W.imag.mean()<0: W=np.conj(W)
    return W
V=[np.array([complex(*p) for p in pts]) for d,pts in res]
D=[d for d,pts in res]
# representatives of the 4 clusters at threshold 0.2
reps=[0,1,6,10]
fig=plt.figure(figsize=(13.0,6.4))
for k,i in enumerate(reps):
    ax=fig.add_subplot(2,2,k+1)
    W=norm_end(V[i])
    ax.plot(W.real,W.imag,'-o',color=['#1f77b4','#d62728','#2ca02c','#9467bd'][k],lw=1.1,ms=3.0)
    ax.set_aspect('equal'); ax.grid(alpha=.18)
    ax.set_title('shape %d: distortion %.5f'%(k+1,D[i]),fontsize=9.5)
    ax.tick_params(labelsize=7.5)
fig.suptitle(r'$m=40$: four macroscopically distinct optima (pairwise sup-distance $>0.2$), '
             r'$v_1\mapsto0$, $v_m\mapsto1$',fontsize=10.5)
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig('../figures/koch-m40-shapes.png',dpi=120,bbox_inches='tight',facecolor='white',pil_kwargs={"optimize":True})
import os
print("bytes:",os.path.getsize('../figures/koch-m40-shapes.png'))
print("\npairwise distances between the four representatives:")
def ip(W,n=1201):
    m=len(W); t=np.linspace(0,1,m); tt=np.linspace(0,1,n)
    return np.interp(tt,t,W.real)+1j*np.interp(tt,t,W.imag)
def variants(W):
    out=[]
    for a in (W,np.conj(W)):
        for rev in (False,True):
            out.append(ip(norm_end(a[::-1] if rev else a)))
    return out
for a in range(4):
    print("   shape %d:"%(a+1), " ".join("%.2f"%min(np.max(np.abs(x-y)) for x in variants(V[reps[a]]) for y in variants(V[reps[b]])) for b in range(4)))
