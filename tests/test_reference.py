import h5py
import numpy as np
import pytest

import radoptlib


def write_sol(path, x, dual, met):
    with h5py.File(path, 'w') as f:
        f.create_dataset('format_note', data='reference solution')
        f.create_dataset('x', data=x)
        f.create_dataset('dual', data=dual)
        f.create_dataset('solver', data='test solver')
        for k, v in met.items():
            f.create_dataset(k, data=float(v))


def test_load_solution_by_path_and_name(tmp_path, monkeypatch, qp_file, qp):
    p = radoptlib.load(qp_file)
    met = radoptlib.verify_solution(p, qp['x'], qp['y'])
    path = tmp_path / 'rt-qp-9999-sol-warmip.h5'
    write_sol(path, qp['x'], qp['y'], met)

    sol = radoptlib.load_solution(path)
    assert sol['name'] == 'rt-qp-9999' and sol['solver'] == 'test solver'
    assert np.array_equal(sol['x'], qp['x']) and np.array_equal(sol['dual'], qp['y'])
    assert sol['objective'] == met['objective'] and sol['gap'] == met['gap']
    assert radoptlib.verify_solution(p, sol['x'], sol['dual'], tol=1e-9)['ok']

    monkeypatch.setenv('RADOPTLIB_DATA', str(tmp_path))
    assert np.array_equal(radoptlib.load_solution('rt-qp-9999')['x'], qp['x'])            # default solver warmip
    write_sol(tmp_path / 'rt-qp-9999-sol-mosek.h5', qp['x'] + 1, qp['y'], met)
    assert np.array_equal(radoptlib.load_solution('rt-qp-9999', solver='mosek')['x'], qp['x'] + 1)


def test_load_solution_missing(tmp_path, monkeypatch):
    monkeypatch.setenv('RADOPTLIB_DATA', str(tmp_path))
    with pytest.raises(FileNotFoundError):
        radoptlib.load_solution('rt-qp-0001')