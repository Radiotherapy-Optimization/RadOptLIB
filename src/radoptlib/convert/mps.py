"""Write an instance as a free-format QPS file (MPS + QUADOBJ, the Maros-Meszaros convention).

Model written:  minimize 0.5 x'Qx + q'x + c  subject to  Gx <= h,  x free.
  - ROWS: objective row 'obj' (N) and one 'L' row per row of G, named c1..cm
  - COLUMNS: variables x1..xn with their q entry and G column
  - RHS: h, and -c on the objective row (the standard way MPS stores an objective constant)
  - BOUNDS: every variable free (MPS defaults to x >= 0 otherwise)
  - QUADOBJ: the lower triangle of Q; readers apply the 0.5 factor
Most QP solvers (Gurobi, CPLEX, MOSEK, HiGHS, COPT, ...) read this format.
The file is written in streaming blocks; large instances give multi-GB text, hence .mps.gz by default.
"""
import gzip
from pathlib import Path

import numpy as np
import scipy.sparse as sp

_BLOCK = 1_000_000          # entries formatted per write


def _num(v):
    return repr(float(v))


def _write_triplets(out, first, second, vals, fmt):
    for s in range(0, len(vals), _BLOCK):
        e = s + _BLOCK
        out.write(''.join(fmt.format(a, b, _num(v)) for a, b, v in zip(first[s:e], second[s:e], vals[s:e])))


def to_mps(problem, path=None, compress=True):
    """Write `problem` (from radoptlib.load) to QPS. Returns the output path.

    path defaults to '<name>.mps.gz' (or '.mps' with compress=False) in the current folder."""
    if path is None:
        path = '{}.mps{}'.format(problem.name, '.gz' if compress else '')
    path = Path(path)
    opener = (lambda p: gzip.open(p, 'wt', compresslevel=6, encoding='ascii')) if str(path).endswith('.gz') \
        else (lambda p: open(p, 'w', encoding='ascii'))
    n, m = problem.n, problem.m
    # objective row on top of G -> one CSC matrix gives every column's entries in order
    A = sp.vstack([sp.csr_matrix(np.asarray(problem.q, dtype=np.float64).reshape(1, n)), problem.G],
                  format='csc')
    A.eliminate_zeros()
    A.sort_indices()
    cols = np.repeat(np.arange(1, n + 1), np.diff(A.indptr))
    row_names = np.array(['obj'] + ['c{}'.format(i) for i in range(1, m + 1)], dtype=object)
    L = sp.tril(problem.Q, format='csc')
    L.eliminate_zeros()
    L.sort_indices()
    q_cols = np.repeat(np.arange(1, n + 1), np.diff(L.indptr))
    with opener(path) as out:
        out.write('NAME {}\nROWS\n N  obj\n'.format(problem.name))
        for s in range(0, m, _BLOCK):
            out.write(''.join(' L  c{}\n'.format(i) for i in range(s + 1, min(s + _BLOCK, m) + 1)))
        out.write('COLUMNS\n')
        _write_triplets(out, cols, row_names[A.indices], A.data, ' x{}  {}  {}\n')
        out.write('RHS\n')
        if problem.c != 0:
            out.write(' RHS  obj  {}\n'.format(_num(-problem.c)))
        h = np.asarray(problem.h, dtype=np.float64)
        nz = np.flatnonzero(h)
        _write_triplets(out, ['RHS'] * len(nz), ['c{}'.format(i + 1) for i in nz], h[nz], ' {}  {}  {}\n')
        out.write('BOUNDS\n')
        for s in range(0, n, _BLOCK):
            out.write(''.join(' FR BND  x{}\n'.format(j) for j in range(s + 1, min(s + _BLOCK, n) + 1)))
        out.write('QUADOBJ\n')
        _write_triplets(out, q_cols, L.indices + 1, L.data, ' x{}  x{}  {}\n')
        out.write('ENDATA\n')
    return path