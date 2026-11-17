import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

sel = json.load(open('sel9b.json'))
cols = ['#1f77b4','#d62728','#2ca02c','#9467bd','#ff7f0e','#8c564b']
al = np.log(3)/np.log(4)
P = [np.array([complex(p[0],p[1]) for p in pts]) for _, pts in sel]

fig = plt.figure(figsize=(12.4, 5.2))
ax  = fig.add_axes([0.045, 0.115, 0.40, 0.80])
axz = fig.add_axes([0.545, 0.115, 0.255, 0.80])
leg = fig.add_axes([0.83, 0.115, 0.16, 0.80]); leg.axis('off')

for W, c, (sp, _) in zip(P, cols, sel):
    ax.plot(W.real, W.imag, '-', color=c, lw=1.0, alpha=.55)
    ax.plot(W.real, W.imag, 'o', color=c, ms=5.0, mec='white', mew=.5)
ax.set_aspect('equal'); ax.grid(alpha=.22)
ax.set_xlabel('normalised real part'); ax.set_ylabel('normalised imaginary part')
ax.set_title(r'$m=9$, $\alpha=\log_43$: six distinct near-optima'
             '\n' r'each normalised by $v_1\mapsto0,\ v_2\mapsto1,\ \operatorname{Im}v_3>0$',
             fontsize=10)
W0 = P[0]
for k in range(9):
    ax.annotate(str(k+1), (W0[k].real, W0[k].imag), textcoords='offset points',
                xytext=(6, 3), fontsize=9, color='0.2')

axz.add_patch(Circle((0,0), 2**al, fill=False, ls=':', ec='0.55', lw=1))
axz.add_patch(Circle((1,0), 1.0, fill=False, ls=':', ec='0.55', lw=1))
for W, c in zip(P, cols):
    axz.plot(W.real, W.imag, '-', color=c, lw=1.2, alpha=.8)
    axz.plot(W.real[:3], W.imag[:3], 'o', color=c, ms=9, mec='white', mew=.6)
axz.plot([0,1],[0,0],'k-',lw=.8)
axz.set_xlim(-0.45, 2.25); axz.set_ylim(-0.95, 1.05); axz.set_aspect('equal'); axz.grid(alpha=.22)
axz.set_title('zoom: the first two points coincide,\nthe third within $0.33$', fontsize=9.5)
axz.set_xlabel(r'dotted: $|z|=2^{\alpha}$ and $|z-1|=1$', fontsize=8.5)

leg.text(0, 0.99, 'distortion', fontsize=9.5, color='0.2')
y = 0.925
for (sp, _), c in zip(sel, cols):
    leg.text(0, y, '%.5f' % sp, fontsize=9.5, color=c); y -= 0.055
leg.text(0, y-0.05, 'spread of $v_k$ across\nthe six (normalised):', fontsize=8.8, color='0.25')
y -= 0.155
sp_k = [max(abs(P[i][k]-P[j][k]) for i in range(6) for j in range(6)) for k in range(9)]
for k in range(9):
    leg.text(0, y, r'$v_%d$: %.2f' % (k+1, sp_k[k]), fontsize=8.6, color='0.3'); y -= 0.045
fig.savefig('../figures/koch-optima-m9.png', dpi=140, bbox_inches='tight', facecolor='white')
print('saved'); print('spreads:', ' '.join('%.2f'%s for s in sp_k))
