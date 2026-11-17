import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
al = np.log(3)/np.log(4); m = 20
sel = json.load(open('m20.json'))[:5]
cols = ['#1f77b4','#d62728','#2ca02c','#9467bd','#ff7f0e']
P = [np.array([complex(*p) for p in s['pts']]) for s in sel]

fig = plt.figure(figsize=(13.2, 5.6))
ax  = fig.add_axes([0.035, 0.10, 0.34, 0.82])
axb = fig.add_axes([0.415, 0.10, 0.34, 0.82])
axr = fig.add_axes([0.79, 0.115, 0.19, 0.79])
leg = fig.add_axes([0.0, 0.0, 0.001, 0.001]); leg.axis('off')

for W, c, s in zip(P, cols, sel):
    ax.plot(W.real, W.imag, '-', color=c, lw=.9, alpha=.6)
    ax.plot(W.real, W.imag, 'o', color=c, ms=4.2, mec='white', mew=.4)
    i,j = s['rmax']
    ax.plot([W[i].real, W[j].real], [W[i].imag, W[j].imag], '--', color=c, lw=1.5, alpha=.9)
ax.set_aspect('equal'); ax.grid(alpha=.2)
ax.set_title(r'five near-optima, $m=20$, $\alpha=\log_43$'+'\n'+r'normalised $v_1\mapsto0,\ v_2\mapsto1$; dashed = max-ratio pair', fontsize=9.5)
ax.set_xlabel('normalised real part'); ax.set_ylabel('normalised imaginary part')

s = sel[0]; W = P[0]
axb.plot(W.real, W.imag, '-', color='0.55', lw=1.0)
axb.plot(W.real, W.imag, 'o', color='0.25', ms=4.6, zorder=5)
for k in range(m):
    axb.annotate(str(k+1), (W[k].real, W[k].imag), textcoords='offset points',
                 xytext=(4,3), fontsize=7.6, color='0.3')
for key, col, lab in [('rmax','#d62728','max'), ('rmin','#1f77b4','min')]:
    i,j = s[key]
    axb.plot([W[i].real, W[j].real], [W[i].imag, W[j].imag], '--', color=col, lw=2.0, zorder=6,
             label=r'%s ratio $%.4f$ at $(v_{%d},v_{%d})$' % (lab, s[key+'v'], i+1, j+1))
axb.legend(loc='upper right', fontsize=7.8, frameon=False)
axb.set_aspect('equal'); axb.grid(alpha=.2)
axb.set_title('best sample: dashed arcs realising the min and max ratio', fontsize=9.5)

r = (np.abs(W[:,None]-W[None,:])/np.abs(np.arange(m)[:,None]-np.arange(m)[None,:])**al)
for d in range(1,m):
    vals = [r[k,k+d] for k in range(m-d)]
    axr.plot([d]*len(vals), vals, '.', color='0.35', ms=2.6)
axr.axhline(s['rmaxv'], color='#d62728', ls='--', lw=1.2)
axr.axhline(s['rminv'], color='#1f77b4', ls='--', lw=1.2)
axr.set_ylim(0.9, 1.36); axr.set_xlabel(r'separation $|i-j|$', fontsize=8.5)
axr.set_ylabel(r'$r_{ij}$', fontsize=8.5); axr.grid(alpha=.2); axr.tick_params(labelsize=7.5)
axr.set_title('ratios by separation', fontsize=9)
fig.savefig('../figures/koch-optima-m20.png', dpi=140, bbox_inches='tight', facecolor='white')
print('saved')
fin = [ (s['dist'], (s['rmin'][0]+1,s['rmin'][1]+1,s['rminv']), (s['rmax'][0]+1,s['rmax'][1]+1,s['rmaxv'])) for s in sel]
for d,a,b in fin: print("  %.5f   min %s   max %s" % (d,a,b))
