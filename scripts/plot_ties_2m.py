import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
al=np.log(3)/np.log(4)
def load(m,f):
    V=np.load(f); V=V-V[0]; V=V/abs(V[1]-V[0]); V=(V-V[0])/(V[-1]-V[0])
    if V.imag.mean()<0: V=np.conj(V)
    R=np.abs(V[:,None]-V[None,:])/np.abs(np.arange(m)[:,None]-np.arange(m)[None,:])**al
    np.fill_diagonal(R,np.nan)
    return V,R,np.nanmax(R),np.nanmin(R)
def count(R,hi,lo,tol):
    A=[(i,j) for i in range(len(R)) for j in range(i+1,len(R)) if R[i,j]>=hi*(1-tol)]
    B=[(i,j) for i in range(len(R)) for j in range(i+1,len(R)) if R[i,j]<=lo*(1+tol)]
    return A,B
D={}
for m,f in ((40,'m40_smooth.npy'),(200,'m200_cont.npy')):
    V,R,hi,lo=load(m,f)
    target=2*m-3                      # the bound
    a,b=1e-7,1e-1
    for _ in range(60):               # bisect the tolerance so that the count ~ 2m-3
        mid=(a*b)**0.5; A,B=count(R,hi,lo,mid)
        if len(A)+len(B) < target: a=mid
        else: b=mid
    tol=(a*b)**0.5; A,B=count(R,hi,lo,tol)
    D[m]=(V,tol,A,B,hi/lo)
    print("m=%3d  distortion %.6f  tol %.3e -> %d max + %d min = %d ties  (2m-3 = %d)"
          %(m,hi/lo,tol,len(A),len(B),len(A)+len(B),target))
    print("       normalised max gaps:", sorted(round((j-i)/m,3) for i,j in A)[:14])
    print("       normalised min gaps:", sorted(round((j-i)/m,3) for i,j in B)[:14])
col={40:'#1f77b4',200:'#9467bd'}
fig=plt.figure(figsize=(12.6,5.5))
ax=fig.add_axes([0.045,0.10,0.40,0.80]); axg=fig.add_axes([0.545,0.10,0.42,0.80])
for m in (40,200):
    V,tol,A,B,d=D[m]
    ax.plot(V.real,V.imag,'-',color=col[m],lw=1.0,alpha=.85,label='$m=%d$, dist. %.4f'%(m,d))
    for (i,j) in A:
        ax.plot([V[i].real,V[j].real],[V[i].imag,V[j].imag],'-',color='#d62728',
                lw=0.7 if m==40 else 0.4,alpha=.5 if m==40 else .22)
    for (i,j) in B:
        ax.plot([V[i].real,V[j].real],[V[i].imag,V[j].imag],'-',color='#2ca02c',
                lw=0.7 if m==40 else 0.4,alpha=.5 if m==40 else .22)
ax.set_aspect('equal'); ax.grid(alpha=.15); ax.legend(fontsize=8.2,loc='lower right')
ax.set_title('near-ties at the tolerance giving $\\approx 2m-3$ pairs\n'
             '($m=40$: %.0e, %d pairs; $m=200$: %.0e, %d pairs)'
             %(D[40][1],len(D[40][2])+len(D[40][3]),D[200][1],len(D[200][2])+len(D[200][3])),fontsize=9)
ax.set_xlabel('real'); ax.set_ylabel('imaginary')
for m,mk in ((40,'o'),(200,'s')):
    V,tol,A,B,d=D[m]
    axg.plot([(i+1)/m for i,j in A],[(j+1)/m for i,j in A],mk,color='#d62728',ms=4.2,mfc='none' if m==200 else '#d62728',mew=1.1 if m==200 else .8,label='$m=%d$ max'%m,alpha=.75 if m==200 else 1)
    axg.plot([(i+1)/m for i,j in B],[(j+1)/m for i,j in B],mk,color='#2ca02c',ms=4.2,mfc='none' if m==200 else '#2ca02c',mew=1.1 if m==200 else .8,label='$m=%d$ min'%m,alpha=.75 if m==200 else 1)
axg.plot([0,1],[0,1],color='0.85',lw=1)
axg.set_xlim(0,1); axg.set_ylim(0,1); axg.grid(alpha=.15); axg.legend(fontsize=7.6,loc='upper left',ncol=2)
axg.set_xlabel('$i/m$'); axg.set_ylabel('$j/m$')
axg.set_title('the tied pairs in normalised parameter',fontsize=9.5)
fig.savefig('../figures/koch-ties-2m.png',dpi=125,bbox_inches='tight',facecolor='white',pil_kwargs={"optimize":True})
import os; print("bytes:",os.path.getsize('../figures/koch-ties-2m.png'))
