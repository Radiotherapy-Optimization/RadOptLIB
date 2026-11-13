# QP instances (`qp`)

## Formulation

```
minimize    ½ x'Qx + q'x + c
subject to  Gx ≤ h
            l ≤ x ≤ u
```

x ∈ ℝⁿ is continuous, Q (n × n) is symmetric positive semidefinite, so every instance is a convex
quadratic program; G is m × n. Bounds on variables are stored as l and u (entries may be
−∞ / +∞).

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
| `small` | below 100,000 | 25 |
| `medium` | 100,000 to 499,999 | 97 |
| `large` | 500,000 and above | 20 |

```python
radoptlib.select(size='small')                     # one class
radoptlib.download('qp', size=['small', 'medium'])  # several classes
radoptlib.select(max_n_vars=50000, min_n_ineq_cons=100000)         # or exact limits on any numeric column
```

## Instance sizes (release 1.0)

| Size                                | Range |
|-------------------------------------|---|
| n_vars                              | 5,513 – 264,034 |
| n_ineq_cons                         | 25,954 – 421,192 |
| nonzeros of Q                       | 0.33 million – 118 million |
| nonzeros of G                       | 1.3 million – 314 million |
| file size                           | 15 MB – 2.0 GB (82 GB in total) |

## Manifest columns (`manifest/qp.csv`)

| Column | Meaning |
|---|---|
| `name` | instance name, `rt-qp-####` |
| `problem_class` | `qp` |
| `problem_type` | `convex continuous QP` |
| `patient` | anonymous patient label, for example `lung-0001`; the same label for all instances of one patient |
| `size` | `small`, `medium` or `large` (see Size classes) |
| `n_vars` | number of variables n |
| `n_ineq_cons` | number of rows m of G |
| `n_eq_cons` | number of equality constraints (always 0 in this class) |
| `nnz_Q`, `nnz_G` | stored nonzeros of Q and G |
| `file_size_bytes` | size of the h5 file |
| `disease_site` | for example `lung`, `oligomets`, `prostate` (more values may be added) |
| `modality` | treatment modality, for example `vmat` |
| `protocol` | prescription of the plan, `<site>_<dose per fraction>gy_<fractions>fx`, for example `prostate_5gy_5fx` |
| `type` | 0, 1 or 2: the stage of the hierarchical optimization (see Origin) |
| `url` | download address |
| `md5` | md5 checksum of the h5 file |
| `sol_url`, `sol_md5` | download address and md5 checksum of the Warm-IP reference solution file (see Reference solution files) |
| `added_in_version` | release in which the instance first appeared |
| `note` | free text, usually empty |

## File format (HDF5, format version 1.0)

Datasets, in the order they are stored in the file:

| Dataset | Content |
|---|---|
| `Q_data`, `Q_indices`, `Q_indptr`, `Q_shape` | Q in CSR form (float64 values, integer indices, shape (n, n)) |
| `q` (n) | float64 vector |
| `c` | float64 scalar |
| `G_data`, `G_indices`, `G_indptr`, `G_shape` | G in CSR form (shape (m, n)) |
| `h` (m) | float64 vector |
| `l` (n), `u` (n) | float64 bounds on x; `-inf` / `inf` where a variable has no bound |
| `data_note` | text, one `key: value` line per root attribute (below), for viewers that do not show attributes |
| `format_note` | text: the formulation and this layout |

Root attributes: `name`, `problem_class`, `problem_type`, `disease_site`, `modality`, `type`, `patient`,
`protocol` (as in the manifest), `library_version`, `format_version`, `license` and `repo_url`.

Some viewers (for example HDFView) list datasets alphabetically by default; the stored order above is
what h5py's `keys()` returns. Q is stored in full (both triangles). Reading a file without the package:

```python
import h5py, scipy.sparse as sp

with h5py.File('rt-qp-0001.h5', 'r') as f:
    csr = lambda k: sp.csr_matrix((f[k + '_data'][()], f[k + '_indices'][()], f[k + '_indptr'][()]),
                                  shape=tuple(f[k + '_shape'][()]))
    Q, q, c = csr('Q'), f['q'][()], float(f['c'][()])
    G, h, l, u = csr('G'), f['h'][()], f['l'][()], f['u'][()]
```

## Reference solution files

`rt-qp-####-sol-<solver>.h5` (for example `rt-qp-0001-sol-warmip.h5`):

| Dataset | Content |
|---|---|
| `x` (n) | primal solution |
| `dual_ineq_cons` (m) | multipliers of Gx ≤ h (≥ 0) |
| `dual_lower_bound` (n), `dual_upper_bound` (n) | multipliers of x ≥ l and x ≤ u (≥ 0; zero where the bound is infinite) |
| `solver` | solver name and link |
| `format_note` | text describing this layout |

Root attributes: `name` (the instance the solution belongs to), `solver`, `format_version`,
`library_version`, `license` and `repo_url`.

The multipliers satisfy Qx + q + G'·dual_ineq_cons − dual_lower_bound + dual_upper_bound = 0. The files
store no objective value or residuals, because those depend on how they are defined; compute them with
`radoptlib.verify_solution`:

```python
sol = radoptlib.load_solution('rt-qp-0001')
radoptlib.verify_solution(radoptlib.load('rt-qp-0001'), sol['x'], sol['dual_ineq_cons'],
                          sol['dual_lower_bound'], sol['dual_upper_bound'])
```
