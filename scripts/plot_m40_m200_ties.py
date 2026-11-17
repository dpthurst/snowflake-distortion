import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
al=np.log(3)/np.log(4)
def load(m,f):
    V=np.load(f); V=V-V[0]; V=V/abs(V[1]-V[0])
    V=(V-V[0])/(V[-1]-V[0])            # endpoint normalisation
    if V.imag.mean()<0: V=np.conj(V)
    R=np.abs(V[:,None]-V[None,:])/np.abs(np.arange(m)[:,None]-np.arange(m)[None,:])**al
    np.fill_diagonal(R,np.nan)
    return V,R,np.nanmax(R),np.nanmin(R)
tol=1e-5
D={}
for m,f in ((40,'m40_smooth.npy'),(200,'m200_cont.npy')):
    V,R,hi,lo=load(m,f)
    A=[(i,j) for i in range(m) for j in range(i+1,m) if R[i,j]>=hi*(1-tol)]
    B=[(i,j) for i in range(m) for j in range(i+1,m) if R[i,j]<=lo*(1+tol)]
    D[m]=(V,A,B,hi/lo)
    print("m=%d: %d max, %d min ; normalised max gaps %s ; min gaps %s" %
          (m,len(A),len(B),[round((j-i)/m,3) for i,j in A],[round((j-i)/m,3) for i,j in B]))
col={40:'#1f77b4',200:'#9467bd'}
fig=plt.figure(figsize=(12.4,5.4))
ax=fig.add_axes([0.045,0.10,0.40,0.80]); axg=fig.add_axes([0.545,0.10,0.42,0.80])
for m in (40,200):
    V,A,B,d=D[m]
    ax.plot(V.real,V.imag,'-',color=col[m],lw=0.9,alpha=.75,label='$m=%d$ (dist. %.4f)'%(m,d))
    for (i,j) in A:
        ax.plot([V[i].real,V[j].real],[V[i].imag,V[j].imag],'--',color='#d62728',lw=1.5 if m==40 else 0.9,alpha=.9 if m==40 else .55)
    for (i,j) in B:
        ax.plot([V[i].real,V[j].real],[V[i].imag,V[j].imag],'--',color='#2ca02c',lw=1.5 if m==40 else 0.9,alpha=.9 if m==40 else .55)
ax.set_aspect('equal'); ax.grid(alpha=.18); ax.legend(fontsize=8.4,loc='lower right')
ax.set_title(r'$m=40$ (heavy) and $m=200$ (light): chords within $10^{-5}$'+'\nred = max side, green = min side',fontsize=9.5)
ax.set_xlabel('real'); ax.set_ylabel('imaginary')
for m in (40,200):
    V,A,B,d=D[m]
    axg.plot([(i+1)/m for i,j in A],[(j+1)/m for i,j in A],'o',color='#d62728',
             ms=9 if m==40 else 6, mfc='none' if m==200 else '#d62728',mew=1.6,label='$m=%d$ max'%m)
    axg.plot([(i+1)/m for i,j in B],[(j+1)/m for i,j in B],'s',color='#2ca02c',
             ms=9 if m==40 else 6, mfc='none' if m==200 else '#2ca02c',mew=1.6,label='$m=%d$ min'%m)
axg.plot([0,1],[0,1],color='0.85',lw=1)
axg.set_xlim(0,1); axg.set_ylim(0,1); axg.grid(alpha=.18); axg.legend(fontsize=8,loc='upper left')
axg.set_xlabel('$i/m$'); axg.set_ylabel('$j/m$')
axg.set_title('near-tied pairs in normalised parameter\n(filled = $m=40$, open = $m=200$)',fontsize=9.5)
fig.savefig('../figures/koch-m40-m200-ties.png',dpi=125,bbox_inches='tight',facecolor='white',pil_kwargs={"optimize":True})
import os; print("bytes:",os.path.getsize('../figures/koch-m40-m200-ties.png'))
