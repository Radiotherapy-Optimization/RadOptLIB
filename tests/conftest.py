"""Synthetic test instances (no real data): a small convex QP with a known optimal primal/dual pair."""
import h5py
import numpy as np
import pytest
import scipy.sparse as sp

FORMAT_NOTE = "minimize 0.5*x'Qx + q'x + c subject to Gx <= h."


def make_qp(n=12, m=20, seed=0):
    """Q PSD; x*, y* >= 0 chosen first, then q and h built so (x*, y*) satisfies the KKT conditions exactly."""
    rng = np.random.default_rng(seed)
    B = sp.random(n, n, density=0.3, random_state=seed, format='csr')
    Q = (B.T @ B + 0.1 * sp.eye(n)).tocsr()
    G = sp.random(m, n, density=0.3, random_state=seed + 1, format='csr')
    G = (G + sp.eye(m, n)).tocsr()                      # no empty rows
    x = rng.standard_normal(n)
    active = rng.random(m) < 0.4
    y = np.where(active, rng.random(m) + 0.1, 0.0)      # y > 0 only on active rows
    slack = np.where(active, 0.0, rng.random(m) + 0.1)
    h = G @ x + slack
    q = -(Q @ x) - G.T @ y
    return {'Q': Q, 'q': q, 'c': 1.5, 'G': G, 'h': h, 'x': x, 'y': y}


def write_h5(path, d, name='rt-qp-9999'):
    note = '\n'.join(['name: ' + name, 'problem_class: qp', 'disease_site: lung', 'modality: vmat', 'type: 1',
                      'protocol: lung_18gy_3fx', 'repo_url: https://example.org', 'format_version: 1.0'])
    with h5py.File(path, 'w') as f:
        f.create_dataset('format_note', data=FORMAT_NOTE)
        f.create_dataset('data_note', data=note)
        for k in ('Q', 'G'):
            M = d[k].tocsr()
            M.sort_indices()
            f.create_dataset(k + '_data', data=M.data, compression='gzip', shuffle=True)
            f.create_dataset(k + '_indices', data=M.indices, compression='gzip', shuffle=True)
            f.create_dataset(k + '_indptr', data=M.indptr, compression='gzip', shuffle=True)
            f.create_dataset(k + '_shape', data=np.asarray(M.shape, dtype=np.int64))
        f.create_dataset('q', data=d['q'])
        f.create_dataset('h', data=d['h'])
        f.create_dataset('c', data=float(d['c']))


@pytest.fixture
def qp():
    return make_qp()


@pytest.fixture
def qp_file(tmp_path, qp):
    p = tmp_path / 'rt-qp-9999.h5'
    write_h5(p, qp)
    return p