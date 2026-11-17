# Configuration data the analysis scripts read

`data/*.npy` holds the configurations behind the figures and Table 1 of the notes. This
directory holds the **other** node sets and optimisation clouds that the scripts in
`scripts/` need as inputs: the earlier or intermediate versions of a configuration, the
branch variants, the perturbation clouds, and the three-point-symmetric candidates. They are
copies, byte for byte, of the files the scripts used to read out of the swarm working
directories, so every script now runs from the repository alone.

| file | m | D | what it is | was |
|---|---|---|---|---|
| `s_m400.npy` | 400 | 1.542271396 | S branch, the version used by the overlap and restart studies | `.ultramath/smwhdqmht/m1/s_m400.npy` |
| `s_m800_plain.npy` | 800 | 1.591051961 | S branch, the version every analysis of m=800 uses (the polished `data/s_m800.npy` came later) | `.ultramath/smwhdqmht/m1/s_m800_plain.npy` |
| `s_m1600_plain.npy` | 1600 | 1.635530900 | S branch; same value as `data/s_m1600.npy` | `.ultramath/smwhdqmht/m1/s_m1600_plain.npy` |
| `s_m1600_kick.npy` | 1600 | 1.635914839 | S branch from a kicked restart | `.ultramath/smwhdqmht/m1/s_m1600_kick.npy` |
| `arc_m400.npy` | 400 | 1.548778919 | arc branch, the version used by the restart study (the polished one is `data/arc400.npy`) | `.ultramath/smwhdqmht/m1/arc_m400.npy` |
| `arc_m800.npy` | 800 | 1.593953536 | arc branch; same value as `data/arc800.npy` | `.ultramath/smwhdqmht/m1/arc_m800.npy` |
| `arc_m1600.npy` | 1600 | 1.631583965 | arc branch; same value as `data/arc1600.npy` | `.ultramath/smwhdqmht/m1/arc_m1600.npy` |
| `S_m800_smooth.npy` | 800 | — | smoothed S branch, for comparison in the alignment tool | `.ultramath/smwhdqmht/m1/S_m800_smooth.npy` |
| `w_smoothS_t2_200.npy`, `w_smoothS_t2_400.npy` | 200, 400 | 1.497314357, 1.543697915 | point-symmetric S branch, second turning level | `.ultramath/smwhdqmht/m2/` |
| `descended_m200.npy` | 200 | 1.503168529 | the descents at m=200; the same value as `data/arc200.npy` | `.scratch/descended_m200.npy` |
| `cert150.npy` | 150 | 1.483370767 | S branch at m=150, outside the ladder | `.scratch/cert150.npy` |
| `final_m40.npz`, `final_m100.npz` | 40, 100 | clouds | 51 and 81 restarts each (`W` node arrays, `Ds` their distortions); the source of the restart and align studies | `.scratch/final_m40.npz`, `.scratch/final_m100.npz` |
| `v2_100best.npy` | 100 | 1.438738049 | best of a 100-start cloud | `.ultramath/smwgx79ya/m8/v2_100best.npy` |
| `coords_3piece_sym.npz` | — | — | parameters of the symmetric three-piece construction | `.ultramath/smwh6zmc9/m4/coords_3piece_sym.npz` |
| `cand_m400.npy`, `cand_m800.npy`, `cand_m1600.npy` | 400, 800, 1600 | — | three-piece candidates, for the alignment tool | `.ultramath/smwh6zmc9/m8/` |

Two consequences worth knowing:

* Where a value here coincides with one in `data/` the two files are still not identical
  arrays: `data/s_m1600.npy` and `s_m1600_plain.npy` agree in $D$ to twelve digits but differ
  by a similarity, and `data/s_m800.npy` (1.590026190) is a *later, better* m=800
  configuration than `s_m800_plain.npy` (1.591051961). The analyses that quote
  `results/shortening-*.md` were run on the files in this directory, which is why they are
  kept rather than merged.
* The reports in `results/` written before this directory existed cite the old
  `.ultramath/...` and `.scratch/...` paths in their prose. The mapping above is the
  translation; the scripts and the recorded `source` fields have been updated to match.

The two shared modules the scripts import also used to live outside the repository;
`scripts/core2.py`, `scripts/slp2.py` and `scripts/bigm.py` are now the copies they import.
