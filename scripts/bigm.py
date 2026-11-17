import numpy as np, time
from scipy.optimize import minimize
al = np.log(3)/np.log(4)
CACHE = {}
def data(m):
    if m not in CACHE:
        idx = np.arange(m)
        CACHE[m] = np.abs(idx[:,None]-idx[None,:]).astype(float)**al
    return CACHE[m]
def V_of(x, m):
    return np.concatenate([[0.0+0j, 1.0+0j], x[:m-2] + 1j*x[m-2:]])
def iu(m): return np.triu_indices(m,1)
def true_dist(V):
    m=len(V); R=(np.abs(V[:,None]-V[None,:])/data(m))[iu(m)]
    return R.max()/R.min(), R
def lse(z):
    mx=z.max(); return mx+np.log(np.exp(z-mx).sum())
def phi_grad(x, m, beta):
    V=V_of(x,m); Dm=data(m)
    R=np.abs(V[:,None]-V[None,:]); np.fill_diagonal(R,1.0)
    with np.errstate(divide='ignore'):
        L=np.log(R)-np.log(Dm)
    idx=iu(m); Lv=L[idx]
    lp=lse(beta*Lv); lm=lse(-beta*Lv)
    val=(lp+lm)/beta
    w=np.exp(beta*Lv-lp); u=np.exp(-beta*Lv-lm)
    c=np.zeros((m,m)); 
    Riu=R[idx]
    c[idx]=(w-u)/Riu**2
    c=c+c.T
    G=np.zeros(m,complex)
    # g_i = sum_j c_ij (v_i - v_j)
    G = V*c.sum(axis=1) - (c@V)
    # note: c@V sums_j c_ij V_j ; and c.sum(axis=1) times V_i
    g=np.concatenate([G[2:].real, G[2:].imag])
    return val, g
def run(m, V0, betas, maxiter=600, verbose=True):
    x=np.concatenate([V0[2:].real,V0[2:].imag])
    best=(true_dist(V_of(x,m))[0], x.copy())
    for beta in betas:
        res=minimize(phi_grad, x, args=(m,beta), jac=True, method='L-BFGS-B',
                     options=dict(maxiter=maxiter, ftol=1e-16, gtol=1e-14))
        x=res.x
        d,r=true_dist(V_of(x,m))
        nact=int((r>=r.max()*(1-1e-3)).sum()+(r<=r.min()*(1+1e-3)).sum())
        if d<best[0]: best=(d,x.copy())
        if verbose: print("      beta %7.0f: %.8f  (ties %d)"%(beta,d,nact), flush=True)
    return best
# validate the gradient
xm=8; IU,_=None,None
IUx=iu(8); rng=np.random.default_rng(0)
x=rng.normal(size=2*8-4)*0.3
val,g=phi_grad(x,8,7.0)
num=np.zeros_like(g)
for k in range(len(x)):
    y=x.copy(); y[k]+=1e-7; num[k]=(phi_grad(y,8,7.0)[0]-val)/1e-7
print("gradient check: max rel err %.2e" % (np.max(np.abs(g-num))/max(1e-12,np.max(np.abs(num)))))
