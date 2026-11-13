# QP instances (`qp`)

The class `qp` contains convex continuous quadratic programs. Other kinds of quadratic program (for
example nonconvex or mixed-integer) will be separate classes with their own names.

- **The problems:** [Formulation](#formulation) · [Origin](#origin) · [Problem sizes](#problem-sizes)
- **The files:** [Instance files](#instance-files) · [Reference solution files](#reference-solution-files)
- **Using the files:** [Reading a file](#reading-a-file) · [Evaluating a solution](#evaluating-a-solution) ·
  [Manifest columns](#manifest-columns-manifestqpcsv)

## Formulation

```
minimize    ½ x'Qx + q'x + c
subject to  Gx ≤ h
            l ≤ x ≤ u
```

x ∈ ℝⁿ is continuous, Q (n × n) is symmetric positive semidefinite and G is m × n. The bounds on the
variables are l and u, whose entries may be −∞ or +∞. Q is not positive definite: in every instance
many variables have no quadratic term.

## Origin

The instances are subproblems from the treatment planning of volumetric modulated arc therapy (VMAT)
plans of real patients, planned with a hierarchical (prioritized) optimization: the clinical goals
are optimized in successive stages, and each stage keeps what the earlier stages achieved while
improving the next goals. The **type** of an instance (0, 1 or 2) is the stage it comes from.
Each patient contributes two or three instances; they share the same anatomy and beam geometry but
differ in the objective terms and constraints of their stage. What the variables and constraints represent is explained under
"What is radiotherapy optimization?" in the [README](../README.md).

## Problem sizes

| Disease site | Patients | Instances | Variables (n) | Constraints (m) | Nonzeros of Q | Nonzeros of G |
|---|---|---|---|---|---|---|
| Lung | 11 | 28 | 6,557 – 131,279 | 39,987 – 172,908 | 0.65 – 16.4 million | 2.9 – 28.0 million |
| Oligometastases | 18 | 54 | 5,513 – 262,463 | 25,954 – 421,192 | 0.33 – 35.3 million | 1.3 – 64.2 million |
| Prostate | 20 | 60 | 77,577 – 264,034 | 177,648 – 412,450 | 26.9 – 118.2 million | 96.6 – 314.2 million |
| **All** | **49** | **142** | 5,513 – 264,034 | 25,954 – 421,192 | 0.33 – 118.2 million | 1.3 – 314.2 million |

The instance files are 15 MB to 2.0 GB each (82 GB in total); the reference solution files are 0.3 to
5.5 MB each (0.3 GB in total).

Every instance has a size class (manifest column `size`) from its total dimension n + m:

| Size | n + m | Instances in 1.0 |
|---|---|---|
| `small` | below 100,000 | 25 |
| `medium` | 100,000 to 499,999 | 97 |
| `large` | 500,000 and above | 20 |

## Instance files

One HDF5 file per instance, `rt-qp-####.h5` (format version 1.0). Datasets, in the order they are
stored in the file:

| Dataset | Content |
|---|---|
| `Q_data`, `Q_indices`, `Q_indptr`, `Q_shape` | Q in CSR form: float64 values, column indices, row pointers (length n + 1) and shape (n, n) |
| `q` (n) | float64 vector |
| `c` | float64 scalar |
| `G_data`, `G_indices`, `G_indptr`, `G_shape` | G in CSR form, shape (m, n) |
| `h` (m) | float64 vector |
| `l` (n), `u` (n) | float64 bounds on x; `-inf` / `inf` where a variable has no bound |
| `data_note` | text, one `key: value` line per root attribute (below), for viewers that do not show attributes |
| `format_note` | text: the formulation and this layout |

Root attributes: `name`, `problem_class`, `problem_type`, `disease_site`, `modality`, `type`, `patient`,
`protocol` (as in the manifest), `library_version`, `format_version`, `license` and `repo_url`.

Column indices and row pointers are 32-bit integers that start at 0, and the column indices are sorted
within each row. Q is stored in full (both triangles).

The files can be looked at without any programming in a free online HDF5 viewer such as
[myHDF5](https://myhdf5.hdfgroup.org): open the page and drop a file onto it. The smallest instance,
`rt-qp-0062.h5` (15 MB), and its reference solution are included in this repository under
[examples/](../examples) as a sample, and these links open them directly in the viewer:
[instance](https://myhdf5.hdfgroup.org/view?url=https%3A%2F%2Fraw.githubusercontent.com%2FRadiotherapy-Optimization%2FRadOptLIB%2Fmain%2Fexamples%2Frt-qp-0062.h5) ·
[reference solution](https://myhdf5.hdfgroup.org/view?url=https%3A%2F%2Fraw.githubusercontent.com%2FRadiotherapy-Optimization%2FRadOptLIB%2Fmain%2Fexamples%2Frt-qp-0062-sol-warmip.h5).

## Reference solution files

One HDF5 file per instance and solver, `rt-qp-####-sol-<solver>.h5` (for example
`rt-qp-0001-sol-warmip.h5`):

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
store no objective value or residuals; they are computed from the vectors as described under
[Evaluating a solution](#evaluating-a-solution).

## Reading a file

The files need no RadOptLIB code: each matrix is built from its stored CSR arrays. The largest
instances (up to 314 million nonzeros in G) need a machine with tens of GB of memory.

<details>
<summary><b>Python</b> (<code>h5py</code> and <code>scipy</code>)</summary>
<br>

```python
import h5py, scipy.sparse as sp

with h5py.File('rt-qp-0001.h5', 'r') as f:
    csr = lambda k: sp.csr_matrix((f[k + '_data'][()], f[k + '_indices'][()], f[k + '_indptr'][()]),
                                  shape=tuple(f[k + '_shape'][()]))
    Q, q, c = csr('Q'), f['q'][()], float(f['c'][()])
    G, h, l, u = csr('G'), f['h'][()], f['l'][()], f['u'][()]

with h5py.File('rt-qp-0001-sol-warmip.h5', 'r') as f:      # a reference solution
    x = f['x'][()]                                         # also dual_ineq_cons, dual_lower_bound, dual_upper_bound
```

</details>

<details>
<summary><b>Julia</b> (packages <code>HDF5</code> and <code>SparseArrays</code>)</summary>
<br>

```julia
using HDF5, SparseArrays

function read_csr(f, key)
    data    = read(f, key * "_data")
    indices = read(f, key * "_indices") .+ 1     # column numbers, 0-based in the file
    indptr  = read(f, key * "_indptr") .+ 1
    m, n    = read(f, key * "_shape")
    At = SparseMatrixCSC(n, m, indptr, indices, data)   # the CSR arrays of A are the CSC arrays of A'
    return sparse(transpose(At))
end

Q, q, c, G, h, l, u = h5open("rt-qp-0001.h5", "r") do f
    read_csr(f, "Q"), read(f, "q"), read(f, "c"), read_csr(f, "G"), read(f, "h"), read(f, "l"), read(f, "u")
end

x = h5read("rt-qp-0001-sol-warmip.h5", "x")      # a reference solution; also dual_ineq_cons,
                                                 # dual_lower_bound, dual_upper_bound
```

</details>

<details>
<summary><b>MATLAB</b> (no toolbox needed)</summary>
<br>

Save the function as `read_csr.m`, or put it at the end of a script.

```matlab
file = 'rt-qp-0001.h5';
Q = read_csr(file, 'Q');
q = h5read(file, '/q');
c = h5read(file, '/c');
G = read_csr(file, 'G');
h = h5read(file, '/h');
l = h5read(file, '/l');     % -Inf where a variable has no lower bound
u = h5read(file, '/u');     % Inf where a variable has no upper bound

x = h5read('rt-qp-0001-sol-warmip.h5', '/x');   % a reference solution; also /dual_ineq_cons,
                                                % /dual_lower_bound, /dual_upper_bound

function A = read_csr(file, key)
    data    = h5read(file, ['/' key '_data']);
    indices = double(h5read(file, ['/' key '_indices'])) + 1;   % column numbers, 0-based in the file
    indptr  = double(h5read(file, ['/' key '_indptr']));
    shape   = double(h5read(file, ['/' key '_shape']));
    rows    = repelem((1:shape(1))', diff(indptr));             % row number of every stored entry
    A = sparse(rows, indices, data, shape(1), shape(2));
end
```

</details>

## Evaluating a solution

`radoptlib.evaluate_solution` reports the objective value and relative residuals of a solution in the
infinity norm, so results from different solvers can be compared on the same scale. With y, z_l and
z_u the multipliers (all ≥ 0) of Gx ≤ h, x ≥ l and x ≤ u:

- `primal_res` = max(‖max(Gx − h, 0)‖, ‖max(l − x, 0)‖, ‖max(x − u, 0)‖) / (1 + max(‖Gx‖, ‖h‖, ‖x‖, ‖l‖, ‖u‖)), with finite bounds only
- `dual_res` = ‖Qx + q + G'y − z_l + z_u‖ / (1 + max(‖Qx‖, ‖q‖, ‖G'y‖, ‖z_l‖, ‖z_u‖))
- `gap` = |pobj − dobj| / (1 + max(|pobj|, |dobj|)), with pobj = ½x'Qx + q'x + c and dobj = −½x'Qx − h'y + l'z_l − u'z_u + c
- `dual_infeas` = how far the multipliers are from being nonnegative, and from zero where a bound is missing

## Manifest columns (`manifest/qp.csv`)

| Column | Meaning |
|---|---|
| `name` | instance name, `rt-qp-####` |
| `problem_class` | `qp` |
| `problem_type` | `convex continuous QP` |
| `patient` | anonymous patient label, for example `lung-0001`; the same label for all instances of one patient |
| `size` | `small`, `medium` or `large` (see [Problem sizes](#problem-sizes)) |
| `n_vars` | number of variables n |
| `n_ineq_cons` | number of rows m of G |
| `n_eq_cons` | number of equality constraints (always 0 in this class) |
| `nnz_Q`, `nnz_G` | stored nonzeros of Q and G |
| `file_size_bytes` | size of the h5 file |
| `disease_site` | for example `lung`, `oligomets`, `prostate` (more values may be added) |
| `modality` | treatment modality, for example `vmat` |
| `protocol` | prescription of the plan, `<site>_<dose per fraction>gy_<fractions>fx`, for example `prostate_5gy_5fx` |
| `type` | 0, 1 or 2: the stage of the hierarchical optimization (see [Origin](#origin)) |
| `url` | download address |
| `md5` | md5 checksum of the h5 file |
| `sol_url`, `sol_md5` | download address and md5 checksum of the Warm-IP reference solution file |
| `added_in_version` | release in which the instance first appeared |
| `note` | free text, usually empty |
