import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
al=np.log(3)/np.log(4)
def kochF(x,depth=40):
    w=np.exp(1j*np.pi/3); d=[1+0j,w,np.conj(w),1+0j]; A=np.array([dk/3 for dk in d]); c=[0+0j]
    for k in range(3): c.append(c[-1]+d[k])
    B=np.array(c)/3
    x=np.atleast_1d(np.asarray(x,float)); n=x.size; y=x.copy(); dig=np.zeros((depth,n),int)
    for i in range(depth):
        t=y*4; k=np.minimum(t.astype(int),3); dig[i]=k; y=t-k
    z=np.zeros(n,complex)
    for i in range(depth-1,-1,-1): z=A[dig[i]]*z+B[dig[i]]
    return z
def norm_end(V):
    W=(V-V[0])/(V[-1]-V[0])
    if W.imag.mean()<0: W=np.conj(W)
    return W
C={}
for m,f in {8:'m8_smooth.npy',20:'m20_smooth.npy',30:'m30_smooth.npy',40:'m40_smooth.npy'}.items():
    V=np.load(f); V=V-V[0]; V=V/abs(V[1]-V[0]); C[m]=norm_end(V)
M=40; W=C[M]; m=M
R=np.abs(W[:,None]-W[None,:])/np.abs(np.arange(m)[:,None]-np.arange(m)[None,:])**al
np.fill_diagonal(R,np.nan)
hi,lo=np.nanmax(R),np.nanmin(R)
A=[(i,j) for i in range(m) for j in range(i+1,m) if R[i,j]>=hi*(1-1e-9)]
B=[(i,j) for i in range(m) for j in range(i+1,m) if R[i,j]<=lo*(1+1e-9)]
print("m=40: ratio in [%.6f, %.6f], distortion %.6f" % (lo,hi,hi/lo))
print("  exact max pairs:", [(i+1,j+1,round(R[i,j],6)) for i,j in A])
print("  exact min pairs:", [(i+1,j+1,round(R[i,j],6)) for i,j in B])

fig=plt.figure(figsize=(12.2,4.7))
ax=fig.add_axes([0.045,0.11,0.37,0.78]); axp=fig.add_axes([0.50,0.11,0.47,0.78])
cols={8:'#1f77b4',20:'#d62728',30:'#2ca02c',40:'#9467bd'}
for mm in sorted(C):
    WW=C[mm]; t=np.linspace(0,1,mm)
    ax.plot(WW.real,WW.imag,'-o',color=cols[mm],lw=1.1,ms=3.2,label='$m=%d$'%mm,alpha=.85)
    axp.plot(t,WW.real,'-',color=cols[mm],lw=1.1); axp.plot(t,WW.imag,'--',color=cols[mm],lw=1.1)
for (i,j) in A:
    ax.plot([W[i].real,W[j].real],[W[i].imag,W[j].imag],'--',color='#d62728',lw=2.2,zorder=8,
            label=r'max chord $(v_{%d},v_{%d})$, gap $%d$, ratio $%.3f\times$min'%(i+1,j+1,j-i,R[i,j]/lo))
for (i,j) in B:
    ax.plot([W[i].real,W[j].real],[W[i].imag,W[j].imag],'--',color='#1f77b4',lw=2.2,zorder=8,
            label=r'min chord $(v_{%d},v_{%d})$, gap $%d$'%(i+1,j+1,j-i))
ax.set_aspect('equal'); ax.grid(alpha=.18); ax.legend(fontsize=7.6,loc='lower right')
ax.set_title('$v_1\\mapsto0$, $v_m\\mapsto1$; extremal chords drawn for $m=40$',fontsize=9.5)
ax.set_xlabel('real'); ax.set_ylabel('imaginary')
axp.grid(alpha=.18); axp.axhline(0,color='0.75',lw=.7)
axp.set_title('coordinates vs $t$ (solid Re, dashed Im)',fontsize=9.5); axp.set_xlabel('$t$')
fig.savefig('../figures/koch-endpoint-overlay.png',dpi=105,bbox_inches='tight',facecolor='white',
            pil_kwargs={"optimize":True})
