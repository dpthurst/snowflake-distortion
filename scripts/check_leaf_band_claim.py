"""Check the claimed constant normalized chord ratio on the three-piece leaf grid."""
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
d = np.load(ROOT / '.ultramath/smwh6zmc9/m4/coords_3piece_sym.npz')
p, z = d['params'], d['coords']
delta = np.diff(p)
keep = delta > 1e-13  # remove only exact shared cell-endpoint parameter records
r = np.abs(np.diff(z)[keep]) / delta[keep] ** (np.log(3)/np.log(4))
print('number of consecutive distinct-parameter pairs:', len(r))
print('rho min:', repr(float(r.min())))
print('rho max:', repr(float(r.max())))
k = int(r.argmax())
indices = np.flatnonzero(keep)
a, b = indices[k], indices[k]+1
print('argmax raw records:', a, b)
print('parameter values:', repr(p[a]), repr(p[b]))
print('parameter gap:', repr(p[b]-p[a]))
print('complex coordinates:', repr(z[a]), repr(z[b]))
assert np.max(np.abs(r)) > 1 + 1e-3
print('The ratio is NOT identically 1 on this saved alpha-similarity arc.')
