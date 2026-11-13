"""Synthetic test instances (no real data): a small convex QP with bounds and a known KKT point."""
import h5py
import numpy as np
import pytest
import scipy.sparse as sp

FORMAT_NOTE = "minimize 0.5*x'Qx + q'x + c subject to Gx <= h and l <= x <= u."
ATTRS = {'name': 'rt-qp-9999', 'problem_class': 'qp', 'problem_type': 'convex continuous QP', 'disease_site': 'lung',
         'modality': 'vmat', 'type': 1, 'patient': 'lung-0001', 'protocol': 'lung_18gy_3fx', 'library_version': '1.0',
         'format_version': '1.0', 'license': 'CC BY 4.0', 'repo_url': 'https://example.org'}


def make_qp(n=12, m=20, seed=0):
    """Q PSD; x*, multipliers >= 0 chosen first, then q, h, l, u built so the KKT conditions hold exactly.

    Variables 0-3: lower bound active (dual_lower_bound > 0); 4-5: upper bound active (dual_upper_bound > 0);
    6-7: both bounds inactive; 8-9: lower bound only (inactive); 10-11: no bounds."""
    rng = np.random.default_rng(seed)
    B = sp.random(n, n, density=0.3, random_state=seed, format='csr')
    Q = (B.T @ B + 0.1 * sp.eye(n)).tocsr()
    Q = ((Q + Q.T) * 0.5).tocsr()
    G = sp.random(m, n, density=0.3, random_state=seed + 1, format='csr')
    G = (G + sp.eye(m, n) + sp.eye(m, n, k=1)).tocsr()   # every row has >= 2 nonzeros
    x = rng.standard_normal(n)
    active = rng.random(m) < 0.4
    y = np.where(active, rng.random(m) + 0.1, 0.0)
    h = G @ x + np.where(active, 0.0, rng.random(m) + 0.1)
    l = np.full(n, -np.inf); u = np.full(n, np.inf)
    zl = np.zeros(n); zu = np.zeros(n)
    l[0:4] = x[0:4]; zl[0:4] = rng.random(4) + 0.1; u[0:4] = x[0:4] + 1.0
    u[4:6] = x[4:6]; zu[4:6] = rng.random(2) + 0.1; l[4:6] = x[4:6] - 2.0
    l[6:8] = x[6:8] - 0.5; u[6:8] = x[6:8] + 0.5
    l[8:10] = x[8:10] - 1.0
    q = -(Q @ x) - G.T @ y + zl - zu
    return {'Q': Q, 'q': q, 'c': 1.5, 'G': G, 'h': h, 'l': l, 'u': u, 'x': x, 'y': y, 'zl': zl, 'zu': zu}


def write_h5(path, d, attrs=ATTRS, bounds=True):
    """bounds=False writes an older-style file (no l/u, no attributes, data_note only)."""
    with h5py.File(path, 'w', track_order=True) as f:
        if bounds:
            for k, v in attrs.items():
                f.attrs[k] = v
        for k in ('Q', 'q', 'c', 'G', 'h', 'l', 'u') if bounds else ('Q', 'q', 'c', 'G', 'h'):
            if k in ('Q', 'G'):
                M = d[k].tocsr()
                M.sort_indices()
                f.create_dataset(k + '_data', data=M.data, compression='gzip', shuffle=True)
                f.create_dataset(k + '_indices', data=M.indices, compression='gzip', shuffle=True)
                f.create_dataset(k + '_indptr', data=M.indptr, compression='gzip', shuffle=True)
                f.create_dataset(k + '_shape', data=np.asarray(M.shape, dtype=np.int64))
            else:
                f.create_dataset(k, data=float(d[k]) if k == 'c' else d[k])
        f.create_dataset('data_note', data='\n'.join('{}: {}'.format(k, v) for k, v in attrs.items()))
        f.create_dataset('format_note', data=FORMAT_NOTE)


@pytest.fixture
def qp():
    return make_qp()


@pytest.fixture
def qp_file(tmp_path, qp):
    p = tmp_path / 'rt-qp-9999.h5'
    write_h5(p, qp)
    return p