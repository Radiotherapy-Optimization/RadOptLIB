import numpy as np
import pytest

import radoptlib
from conftest import write_h5


def test_load_from_path(qp_file, qp):
    p = radoptlib.load(qp_file)
    assert p.name == 'rt-qp-9999'
    assert (p.n, p.m) == qp['G'].shape[::-1]
    assert (p.Q != qp['Q']).nnz == 0 and (p.G != qp['G']).nnz == 0
    assert np.array_equal(p.q, qp['q']) and np.array_equal(p.h, qp['h']) and p.c == 1.5
    assert np.array_equal(p.l, qp['l']) and np.array_equal(p.u, qp['u'])
    assert p.info['disease_site'] == 'lung' and p.info['type'] == 1 and p.info['format_version'] == '1.0'
    x = qp['x']
    assert p.objective(x) == pytest.approx(0.5 * x @ (qp['Q'] @ x) + qp['q'] @ x + 1.5)


def test_load_file_without_bounds(tmp_path, qp):
    path = tmp_path / 'old.h5'
    write_h5(path, qp, bounds=False)
    p = radoptlib.load(path)
    assert np.all(p.l == -np.inf) and np.all(p.u == np.inf)
    assert p.info['type'] == '1'                          # parsed from the data_note text


def test_evaluate_solution_at_optimum(qp_file, qp):
    p = radoptlib.load(qp_file)
    r = radoptlib.evaluate_solution(p, qp['x'], qp['y'], qp['zl'], qp['zu'], tol=1e-9)
    assert r['ok'], r
    assert r['primal_res'] < 1e-12 and r['dual_res'] < 1e-12 and r['gap'] < 1e-12 and r['dual_infeas'] == 0


def test_evaluate_solution_detects_errors(qp_file, qp):
    p = radoptlib.load(qp_file)
    r = radoptlib.evaluate_solution(p, qp['x'] + 1.0, qp['y'], qp['zl'], qp['zu'], tol=1e-6)
    assert not r['ok'] and r['dual_res'] > 1e-6
    r = radoptlib.evaluate_solution(p, qp['x'], -np.ones(p.m), qp['zl'], qp['zu'], tol=1e-6)
    assert not r['ok'] and r['dual_infeas'] == 1.0
    r = radoptlib.evaluate_solution(p, qp['x'], qp['y'], qp['zl'], None, tol=1e-6)    # dual_upper_bound left out
    assert not r['ok'] and r['dual_res'] > 1e-6
    zl = qp['zl'].copy()
    zl[10] = 0.5                                              # multiplier on a missing bound
    assert radoptlib.evaluate_solution(p, qp['x'], qp['y'], zl, qp['zu'])['dual_infeas'] == 0.5


def test_evaluate_solution_bound_violation(qp_file, qp):
    p = radoptlib.load(qp_file)
    x = qp['x'].copy()
    x[0] -= 1e-3                                              # below its lower bound
    assert radoptlib.evaluate_solution(p, x)['primal_res'] > 1e-5
    x = qp['x'].copy()
    x[5] += 1e-3                                              # above its upper bound
    assert radoptlib.evaluate_solution(p, x)['primal_res'] > 1e-5


def test_evaluate_solution_primal_only(qp_file, qp):
    p = radoptlib.load(qp_file)
    r = radoptlib.evaluate_solution(p, qp['x'])
    assert set(r) == {'objective', 'primal_res'} and r['primal_res'] < 1e-12
    bad = radoptlib.evaluate_solution(p, qp['x'] * 0 + 1e6)
    assert bad['primal_res'] > 0


def test_evaluate_solution_wrong_length(qp_file, qp):
    p = radoptlib.load(qp_file)
    with pytest.raises(ValueError):
        radoptlib.evaluate_solution(p, np.zeros(p.n + 1))
    with pytest.raises(ValueError):
        radoptlib.evaluate_solution(p, qp['x'], qp['y'], np.zeros(p.n - 1))