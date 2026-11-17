import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
al = np.log(3)/np.log(4); m = 20
s = json.load(open('m20.json'))[0]
V = np.array([complex(*p) for p in s['pts']])
V = V - V[0]; V = V/abs(V[1]-V[0])
if V[2].imag < 0: V = np.conj(V)
R = np.abs(V[:,None]-V[None,:])/np.abs(np.arange(m)[:,None]-np.arange(m)[None,:])**al
np.fill_diagonal(R, np.nan)
hi = np.nanmax(R); lo = np.nanmin(R)
print("distortion %.5f   max %.6f   min %.6f" % (hi/lo, hi, lo))
for tol in (1e-3, 1e-2):
    A = [(i,j) for i in range(m) for j in range(i+1,m) if R[i,j] >= hi*(1-tol)]
    B = [(i,j) for i in range(m) for j in range(i+1,m) if R[i,j] <= lo*(1+tol)]
    print("tol %.0e : %d near-max pairs, %d near-min pairs" % (tol, len(A), len(B)))

fig = plt.figure(figsize=(13.0, 6.0))
axes = [fig.add_axes([0.04, 0.08, 0.44, 0.80]), fig.add_axes([0.52, 0.08, 0.44, 0.80])]
for ax, tol in zip(axes, (1e-3, 1e-2)):
    ax.plot(V.real, V.imag, '-', color='0.6', lw=0.9, zorder=1)
    ax.plot(V.real, V.imag, 'o', color='0.15', ms=5, zorder=6)
    nA = nB = 0
    for i in range(m):
        for j in range(i+1,m):
            if R[i,j] >= hi*(1-tol):
                ax.plot([V[i].real,V[j].real],[V[i].imag,V[j].imag],'--',color='#d62728',lw=1.6,alpha=.85,zorder=3); nA+=1
            elif R[i,j] <= lo*(1+tol):
                ax.plot([V[i].real,V[j].real],[V[i].imag,V[j].imag],'--',color='#1f77b4',lw=1.6,alpha=.85,zorder=3); nB+=1
    for k in range(m):
        ax.annotate(str(k+1),(V[k].real,V[k].imag),textcoords='offset points',xytext=(5,3),
                    fontsize=7.8,color='0.3')
    ax.set_aspect('equal'); ax.grid(alpha=.2)
    ax.set_title("near-ties within %s:  %d at the max ($\\geq %.4f$),  %d at the min ($\\leq %.4f$)"
                 % (("$0.1\\%$" if tol<5e-3 else "$1\\%$"), nA, hi*(1-tol), nB, lo*(1+tol)), fontsize=9.5)
    ax.set_xlabel('normalised real part')
axes[0].set_ylabel('normalised imaginary part')
fig.suptitle(r"$m=20$, $\alpha=\log_43$, distortion $%.5f$: chords close to the extremes" % (hi/lo), fontsize=11)
fig.savefig('../figures/koch-m20-nearties.png', dpi=140, bbox_inches='tight', facecolor='white')
print("saved")
