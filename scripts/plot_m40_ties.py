import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
al=np.log(3)/np.log(4); m=40
V=np.load('m40_smooth.npy'); V=V-V[0]; V=V/abs(V[1]-V[0])
if V.imag.mean()<0: V=np.conj(V)
R=np.abs(V[:,None]-V[None,:])/np.abs(np.arange(m)[:,None]-np.arange(m)[None,:])**al
np.fill_diagonal(R,np.nan); hi,lo=np.nanmax(R),np.nanmin(R)
tol=1e-5
A=[(i,j,float(R[i,j])) for i in range(m) for j in range(i+1,m) if R[i,j]>=hi*(1-tol)]
B=[(i,j,float(R[i,j])) for i in range(m) for j in range(i+1,m) if R[i,j]<=lo*(1+tol)]
print("m=40, distortion %.8f ; tol %.0e -> %d max-side, %d min-side"%(hi/lo,tol,len(A),len(B)))
print("  max-side:", [(i+1,j+1,round(v/min(v for *_,v in A),6)) for i,j,v in A])
print("  min-side:", [(i+1,j+1,round(v/lo,7)) for i,j,v in B])
fig=plt.figure(figsize=(9.6,6.4))
ax=fig.add_axes([0.07,0.07,0.60,0.86])
ax.plot(V.real,V.imag,'-',color='0.62',lw=1.0,zorder=1)
ax.plot(V.real,V.imag,'o',color='0.15',ms=4.4,zorder=7)
for (i,j,v) in A:
    ax.plot([V[i].real,V[j].real],[V[i].imag,V[j].imag],'--',color='#d62728',lw=1.9,alpha=.9,zorder=4)
for (i,j,v) in B:
    ax.plot([V[i].real,V[j].real],[V[i].imag,V[j].imag],'--',color='#1f77b4',lw=1.9,alpha=.9,zorder=4)
for k in range(m):
    ax.annotate(str(k+1),(V[k].real,V[k].imag),textcoords="offset points",xytext=(4.5,3),
                fontsize=6.6,color='0.25')
ax.set_aspect('equal'); ax.grid(alpha=.18)
ax.set_title(r'$m=40$: chords within $10^{-5}$ of the extreme ratios',fontsize=11)
ax.set_xlabel('real'); ax.set_ylabel('imaginary')
txt=[("max side (ratio $\\geq %.6f$)"%hi, "#d62728")]
txt+= [("  $v_{%d}$–$v_{%d}$"%(i+1,j+1), "#d62728") for i,j,v in A]
txt+= [("min side (ratio $\\leq %.6f$)"%lo, "#1f77b4")]
txt+= [("  $v_{%d}$–$v_{%d}$"%(i+1,j+1), "#1f77b4") for i,j,v in B]
txt+= [("", "0.2"), (r"distortion $%.6f$"%(hi/lo), "0.15"),
       (r"$v_1\mapsto0$, $v_m\mapsto1$", "0.15"),
       (r"gaps: "+", ".join("%d"%(j-i) for i,j,v in A), "#d62728"),
       (r"       "+", ".join("%d"%(j-i) for i,j,v in B), "#1f77b4")]
lx=fig.add_axes([0.685,0.09,0.30,0.80]); lx.axis('off')
y=1.0
for s,c in txt:
    lx.text(0,y,s,fontsize=8.2,color=c,transform=lx.transAxes); y-=0.052
fig.savefig('../figures/koch-m40-ties.png',dpi=125,bbox_inches='tight',facecolor='white',pil_kwargs={"optimize":True})
import os; print("bytes:", os.path.getsize('../figures/koch-m40-ties.png'))
