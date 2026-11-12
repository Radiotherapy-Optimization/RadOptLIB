import numpy as np
import pytest

import radoptlib


def test_load_from_path(qp_file, qp):
    p = radoptlib.load(qp_file)
    assert p.name == 'rt-qp-9999'
    assert (p.n, p.m) == qp['G'].shape[::-1]
    assert (p.Q != qp['Q']).nnz == 0 and (p.G != qp['G']).nnz == 0
    assert np.array_equal(p.q, qp['q']) and np.array_equal(p.h, qp['h']) and p.c == 1.5
    assert p.info['disease_site'] == 'lung' and p.info['type'] == '1' and p.info['format_version'] == '1.0'
    x = qp['x']
    assert p.objective(x) == pytest.approx(0.5 * x @ (qp['Q'] @ x) + qp['q'] @ x + 1.5)


def test_verify_solution_at_optimum(qp_file, qp):
    p = radoptlib.load(qp_file)
    r = radoptlib.verify_solution(p, qp['x'], qp['y'], tol=1e-9)
    assert r['ok'], r
    assert r['primal_res'] < 1e-12 and r['dual_res'] < 1e-12 and r['gap'] < 1e-12 and r['dual_infeas'] == 0


def test_verify_solution_detects_errors(qp_file, qp):
    p = radoptlib.load(qp_file)
    r = radoptlib.verify_solution(p, qp['x'] + 1.0, qp['y'], tol=1e-6)
    assert not r['ok'] and r['dual_res'] > 1e-6
    r = radoptlib.verify_solution(p, qp['x'], -np.ones(p.m), tol=1e-6)
    assert not r['ok'] and r['dual_infeas'] == 1.0


def test_verify_solution_primal_only(qp_file, qp):
    p = radoptlib.load(qp_file)
    r = radoptlib.verify_solution(p, qp['x'])
    assert set(r) == {'objective', 'primal_res'}
    bad = radoptlib.verify_solution(p, qp['x'] * 0 + 1e6)
    assert bad['primal_res'] > 0


def test_verify_solution_wrong_length(qp_file):
    p = radoptlib.load(qp_file)
    with pytest.raises(ValueError):
        radoptlib.verify_solution(p, np.zeros(p.n + 1))