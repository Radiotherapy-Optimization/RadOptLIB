"""to_mps: parse the written QPS with a small independent reader and compare with the original data;
optionally read it with HiGHS (skipped if highspy is not installed)."""
import gzip

import numpy as np
import pytest
import scipy.sparse as sp

import radoptlib


def read_qps(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    sec, rows, cols, A, rhs, bounds, quad, name = None, [], {}, {}, {}, {}, {}, None
    with opener(path, 'rt') as f:
        for line in f:
            t = line.split()
            if not line.startswith(' '):
                sec = t[0]
                if sec == 'NAME':
                    name = t[1]
                continue
            if sec == 'ROWS':
                rows.append((t[0], t[1]))
            elif sec == 'COLUMNS':
                cols.setdefault(t[0], len(cols))
                A[(t[1], t[0])] = float(t[2])
            elif sec == 'RHS':
                rhs[t[1]] = float(t[2])
            elif sec == 'BOUNDS':
                bounds.setdefault(t[2], []).append((t[0], float(t[3]) if len(t) > 3 else None))
            elif sec == 'QUADOBJ':
                quad[(t[0], t[1])] = float(t[2])
    return name, rows, cols, A, rhs, bounds, quad


@pytest.mark.parametrize('compress', [True, False])
def test_mps_roundtrip(tmp_path, qp_file, compress):
    p = radoptlib.load(qp_file)
    path = radoptlib.to_mps(p, tmp_path / ('x.mps.gz' if compress else 'x.mps'))
    name, rows, cols, A, rhs, bounds, quad = read_qps(path)
    n, m = p.n, p.m
    assert name == 'rt-qp-9999'
    assert rows[0] == ('N', 'obj') and len(rows) == m + 1 and all(r[0] == 'L' for r in rows[1:])
    assert len(bounds) == n
    l = np.zeros(n); u = np.full(n, np.inf)                # MPS defaults
    for nm, entries in bounds.items():
        j = int(nm[1:]) - 1
        for kind, v in entries:
            if kind == 'FR':
                l[j], u[j] = -np.inf, np.inf
            elif kind == 'MI':
                l[j] = -np.inf
            elif kind == 'LO':
                l[j] = v
            elif kind == 'UP':
                u[j] = v
            elif kind == 'FX':
                l[j] = u[j] = v
    assert np.array_equal(l, p.l) and np.array_equal(u, p.u)
    q = np.zeros(n); G = np.zeros((m, n)); Q = np.zeros((n, n))
    for (r, c), v in A.items():
        j = int(c[1:]) - 1
        if r == 'obj':
            q[j] = v
        else:
            G[int(r[1:]) - 1, j] = v
    for (a, b), v in quad.items():
        i, j = int(a[1:]) - 1, int(b[1:]) - 1
        assert i <= j                                     # lower triangle: row index >= column index
        Q[i, j] = Q[j, i] = v
    h = np.array([rhs.get('c{}'.format(i + 1), 0.0) for i in range(m)])
    assert np.array_equal(q, p.q)                          # repr(float) round-trips exactly
    assert np.array_equal(G, p.G.toarray())
    assert np.array_equal(Q, p.Q.toarray())
    assert np.array_equal(h, p.h)
    assert -rhs['obj'] == p.c


def test_mps_highs(tmp_path, qp_file, qp):
    highspy = pytest.importorskip('highspy')
    p = radoptlib.load(qp_file)
    path = radoptlib.to_mps(p, tmp_path / 'x.mps', compress=False)
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    assert h.readModel(str(path)) == highspy.HighsStatus.kOk
    h.run()
    x = np.array(h.getSolution().col_value)
    assert h.getInfo().objective_function_value == pytest.approx(p.objective(qp['x']), rel=1e-6)
    assert radoptlib.verify_solution(p, x)['primal_res'] < 1e-7
    assert np.array_equal(np.array(h.getLp().col_lower_), np.maximum(p.l, -highspy.kHighsInf))