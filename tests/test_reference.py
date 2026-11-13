import h5py
import numpy as np
import pytest

import radoptlib


def write_sol(path, qp, shift=0.0):
    with h5py.File(path, 'w', track_order=True) as f:
        f.attrs['name'] = 'rt-qp-9999'
        f.attrs['solver'] = 'test solver'
        f.create_dataset('x', data=qp['x'] + shift)
        f.create_dataset('dual_ineq_cons', data=qp['y'])
        f.create_dataset('dual_lower_bound', data=qp['zl'])
        f.create_dataset('dual_upper_bound', data=qp['zu'])
        f.create_dataset('solver', data='test solver')
        f.create_dataset('format_note', data='reference solution')


def test_load_solution_by_path_and_name(tmp_path, monkeypatch, qp_file, qp):
    p = radoptlib.load(qp_file)
    path = tmp_path / 'rt-qp-9999-sol-warmip.h5'
    write_sol(path, qp)

    sol = radoptlib.load_solution(path)
    assert sol['name'] == 'rt-qp-9999' and sol['solver'] == 'test solver'
    assert set(sol) == {'name', 'solver', 'format_note', 'x', 'dual_ineq_cons', 'dual_lower_bound',
                        'dual_upper_bound'}                                               # no stored residuals
    assert np.array_equal(sol['x'], qp['x']) and np.array_equal(sol['dual_ineq_cons'], qp['y'])
    assert np.array_equal(sol['dual_lower_bound'], qp['zl']) and np.array_equal(sol['dual_upper_bound'], qp['zu'])
    assert radoptlib.verify_solution(p, sol['x'], sol['dual_ineq_cons'], sol['dual_lower_bound'],
                                     sol['dual_upper_bound'], tol=1e-9)['ok']

    monkeypatch.setenv('RADOPTLIB_DATA', str(tmp_path))
    assert np.array_equal(radoptlib.load_solution('rt-qp-9999')['x'], qp['x'])            # default solver warmip
    write_sol(tmp_path / 'rt-qp-9999-sol-mosek.h5', qp, shift=1.0)
    assert np.array_equal(radoptlib.load_solution('rt-qp-9999', solver='mosek')['x'], qp['x'] + 1)


def test_load_solution_missing(tmp_path, monkeypatch):
    monkeypatch.setenv('RADOPTLIB_DATA', str(tmp_path))
    with pytest.raises(FileNotFoundError):
        radoptlib.load_solution('rt-qp-0001')