# QP instances (`qp`)

## Formulation

```
minimize    ½ x'Qx + q'x + c
subject to  Gx ≤ h
```

x ∈ ℝⁿ is continuous, Q (n × n) is symmetric positive semidefinite, so every instance is a convex
quadratic program; G is m × n. There are no other constraints: bounds on variables appear as rows
of G.

## Origin

The instances are subproblems from the treatment planning of volumetric modulated arc therapy (VMAT)
plans of real patients, planned with a hierarchical (prioritized) optimization: the clinical goals
are optimized in successive stages, and each stage keeps what the earlier stages achieved while
improving the next goals. The **type** of an instance (0, 1 or 2) is the stage it comes from.
Instances of the same patient share the same anatomy and beam geometry but differ in the objective
terms and constraints of their stage.

## Size classes

Every instance has a size class (manifest column `size`) from its total dimension n + m. The limits
are fixed, so later instances are classified by the same rule.

| Size | n + m | Instances in 1.0 |
|---|---|---|
| `small` | below 100,000 | 20 |
| `medium` | 100,000 to 499,999 | 91 |
| `large` | 500,000 and above | 31 |

```python
radoptlib.select(size='small')                     # one class
radoptlib.download('qp', size=['small', 'medium'])  # several classes
radoptlib.select(max_n=50000, min_m=100000)         # or exact limits on any numeric column
```

## Instance sizes (release 1.0)

| Size | Range |
|---|---|
| variables n | 5,513 – 264,034 |
| inequality rows m | 34,179 – 726,461 |
| nonzeros of Q | 0.33 million – 118 million |
| nonzeros of G | 1.3 million – 314 million |
| file size | 15 MB – 2.0 GB (82 GB in total) |

## Manifest columns (`manifest/qp.csv`)

| Column | Meaning |
|---|---|
| `name` | instance name, `rt-qp-####` |
| `problem_class` | `qp` |
| `problem_type` | `convex continuous QP` |
| `disease_site` | for example `lung`, `oligomets`, `prostate` (more values may be added) |
| `modality` | treatment modality, for example `vmat` |
| `type` | 0, 1 or 2: the stage of the hierarchical optimization (see Origin) |
| `patient` | anonymous patient label, for example `lung-0001`; the same label for all instances of one patient |
| `protocol` | prescription of the plan, `<site>_<dose per fraction>gy_<fractions>fx`, for example `prostate_5gy_5fx` |
| `size` | `small`, `medium` or `large` (see Size classes) |
| `n`, `m` | number of variables and of rows of G |
| `n_equalities` | number of equality constraints (always 0 in this class) |
| `nnz_Q`, `nnz_G` | stored nonzeros of Q and G |
| `file_size_bytes` | size of the h5 file |
| `url` | download address |
| `added_in_version` | release in which the instance first appeared |
| `note` | free text, usually empty |
| `md5` | md5 checksum of the h5 file |

## File format (HDF5, format version 1.0)

| Dataset | Content |
|---|---|
| `Q_data`, `Q_indices`, `Q_indptr`, `Q_shape` | Q in CSR form (float64 values, integer indices, shape (n, n)) |
| `G_data`, `G_indices`, `G_indptr`, `G_shape` | G in CSR form (shape (m, n)) |
| `q` (n), `h` (m) | float64 vectors |
| `c` | float64 scalar |
| `format_note` | text: the formulation and this layout |
| `data_note` | text, `key: value` lines: name, problem_class, disease_site, modality, type, protocol, repo_url, format_version |

Q is stored in full (both triangles). Reading a file without the package:

```python
import h5py, scipy.sparse as sp

with h5py.File('rt-qp-0001.h5', 'r') as f:
    csr = lambda k: sp.csr_matrix((f[k + '_data'][()], f[k + '_indices'][()], f[k + '_indptr'][()]),
                                  shape=tuple(f[k + '_shape'][()]))
    Q, G, q, h, c = csr('Q'), csr('G'), f['q'][()], f['h'][()], float(f['c'][()])
```

## Reference solution files

`rt-qp-####-sol-<solver>.h5` (for example `rt-qp-0001-sol-warmip.h5`):

| Dataset | Content |
|---|---|
| `x` (n) | primal solution |
| `dual` (m) | multipliers of Gx ≤ h (≥ 0) |
| `objective`, `primal_res`, `dual_res`, `gap`, `dual_infeas` | values computed with `radoptlib.verify_solution` |
| `solver` | solver name and link |
| `format_note` | text describing this layout |
