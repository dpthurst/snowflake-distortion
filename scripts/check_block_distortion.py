"""Exhaustively compare finite distortion of every long contiguous block.

For each one-row ratio R[i,j]=|v_i-v_j|/(j-i)^alpha, first take prefix
extrema over j, then use a one-dimensional sparse range table over i.
This queries all pairs inside a block without forming it. Random literal
submatrix checks validate both the maximum and minimum before the full scan.
The whole interval is included in the count, unlike some earlier reports.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/s_m800.npy'
v = np.load(SOURCE)
assert v.shape == (800,) and np.isfinite(v).all()
N = len(v)
alpha = np.log(3)/np.log(4)
I, J = np.triu_indices(N, 1)
r = np.abs(v[I]-v[J])/(J-I)**alpha
D0 = r.max()/r.min()
R = np.full((N,N),np.nan)
R[I,J] = r


def sparse_rows(a, op):
    # At (i,j), first cover j' in [i+1,j], then the row RMQ below covers i.
    tables = [op.accumulate(a, axis=1)]
    for k in range(1, N.bit_length()):
        h = 1 << (k-1)
        prev = tables[-1]
        new = prev.copy()
        new[:N-2*h+1] = op(prev[:N-2*h+1,:], prev[h:N-h+1,:])
        tables.append(new)
    return tables


hi = sparse_rows(np.where(np.isfinite(R), R, -np.inf), np.maximum)
lo = sparse_rows(np.where(np.isfinite(R), R, np.inf), np.minimum)


def query(tables, L):
    """Return extrema for all blocks [a,a+L], 0 <= a < N-L."""
    k = L.bit_length()-1
    s = 1 << k
    n = N-L
    starts = np.arange(n)
    second = starts+L-s
    op = np.maximum if tables is hi else np.minimum
    return op(tables[k][starts,starts+L], tables[k][second,starts+L])


rng = np.random.default_rng(0)
for _ in range(40):
    a = int(rng.integers(N-2)); b = int(rng.integers(a+1,N))
    block = R[a:b+1,a:b+1]
    L = b-a
    np.testing.assert_allclose([query(hi,L)[a],query(lo,L)[a]],
                               [np.nanmax(block),np.nanmin(block)],rtol=0,atol=0)
print('40 literal RMQ spot checks: PASS')
rows = []
for L in range(160,800):
    D = query(hi,L)/query(lo,L)
    rows.append((L,int(D.argmin()),float(D.min()),float(np.max(np.abs(D-D0))),
                 int(np.sum(np.abs(D-D0)<=D0*1e-9))))
rows = np.array(rows)
for row in rows:
    assert row[3] >= 0
count = int((N-np.arange(160,N)).sum())
assert count == 205120
assert np.isclose(rows[:,3].max(), 4.5523760454635465e-06)
assert int(rows[:,4].sum()) == 50616
result = {
    'source': str(SOURCE.relative_to(ROOT)),
    'interval_count': count,
    'minimum_index_gap': 160,
    'global_D': float(D0),
    'max_absolute_deviation': float(rows[:,3].max()),
    'max_relative_deviation': float((rows[:,3]/D0).max()),
    'relative_tie_tolerance': 1e-9,
    'near_tie_count': int(rows[:,4].sum()),
    'validation': '40 literal random-submatrix max/min checks passed',
    'rows_L_worst_a_min_D_max_abs_error_near_ties': rows.tolist(),
}
Path(ROOT/'results/block-distortion-exhaustive.json').write_text(
    json.dumps(result, indent=2)+'\n'
)
for key in ('interval_count','global_D','max_absolute_deviation',
            'max_relative_deviation','near_tie_count','validation'):
    print(f'{key}: {result[key]}')
