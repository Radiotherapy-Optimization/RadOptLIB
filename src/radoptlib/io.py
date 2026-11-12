"""Manifest lookup, downloads and loading of RadOptLIB instances."""
import csv
import hashlib
import io
import os
import shutil
import urllib.request
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

PROBLEM_CLASSES = ('qp',)
MANIFEST_URL = 'https://raw.githubusercontent.com/Radiotherapy-Optimization/RadOptLIB/main/manifest/{}.csv'
INT_COLUMNS = ('type', 'n', 'm', 'n_equalities', 'nnz_Q', 'nnz_G', 'file_size_bytes')

_manifest_cache = {}


def data_dir():
    """Folder for downloaded files: $RADOPTLIB_DATA if set, else ~/.radoptlib."""
    d = Path(os.environ.get('RADOPTLIB_DATA', Path.home() / '.radoptlib'))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _parse_manifest(text):
    rows = list(csv.DictReader(io.StringIO(text)))
    for r in rows:
        for k in INT_COLUMNS:
            if r.get(k, '') != '':
                r[k] = int(r[k])
    return rows


def _load_manifest(problem_class):
    """Manifest rows of one problem class.

    Order: the repo's own manifest/ folder when running from a checkout, else the copy on GitHub
    (cached in data_dir()/manifest/), else the cached copy when offline."""
    if problem_class in _manifest_cache:
        return _manifest_cache[problem_class]
    local = Path(__file__).resolve().parents[2] / 'manifest' / '{}.csv'.format(problem_class)
    if local.exists():
        text = local.read_text(encoding='utf-8')
    else:
        cached = data_dir() / 'manifest' / '{}.csv'.format(problem_class)
        try:
            with urllib.request.urlopen(MANIFEST_URL.format(problem_class), timeout=30) as resp:
                text = resp.read().decode('utf-8')
            cached.parent.mkdir(parents=True, exist_ok=True)
            cached.write_text(text, encoding='utf-8')
        except OSError:
            if not cached.exists():
                raise RuntimeError('cannot download the {} manifest and no cached copy exists'.format(problem_class))
            text = cached.read_text(encoding='utf-8')
    _manifest_cache[problem_class] = _parse_manifest(text)
    return _manifest_cache[problem_class]


def list_instances(problem_class=None):
    """All manifest rows (dicts), optionally of one problem class."""
    classes = PROBLEM_CLASSES if problem_class is None else (problem_class.lower(),)
    return [r for c in classes for r in _load_manifest(c)]


def info(name):
    """Manifest row of one instance, e.g. info('rt-qp-0001')."""
    for r in list_instances():
        if r['name'] == name:
            return dict(r)
    raise KeyError('unknown instance: {}'.format(name))


def _matches(row, key, value):
    if key not in row and key[:4] in ('min_', 'max_') and key[4:] in row:      # range filter on a number column
        col = key[4:]
        if col not in INT_COLUMNS:
            raise KeyError('range filters need a numeric column, not {}'.format(col))
        return row[col] >= value if key.startswith('min_') else row[col] <= value
    if key not in row:
        raise KeyError('unknown manifest column: {}'.format(key))
    values = value if isinstance(value, (list, tuple, set)) else [value]
    return any(str(row[key]).lower() == str(v).lower() for v in values)


def select(**filters):
    """Names of the instances matching every filter, e.g. select(problem_class='qp', disease_site='lung', type=1).

    Any manifest column can be a filter; a list value means "any of", e.g. size=['small', 'medium'].
    Numeric columns also take inclusive ranges with min_/max_ prefixes, e.g.
    select(max_n=50000, min_m=100000) or select(max_file_size_bytes=500e6)."""
    return [r['name'] for r in list_instances()
            if all(_matches(r, k, v) for k, v in filters.items())]


def _md5(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def _is_cached(path, row):
    if not path.exists() or path.stat().st_size != row['file_size_bytes']:
        return False
    return not row.get('md5') or _md5(path) == row['md5']


def download(name_or_class, disease_site=None, **filters):
    """Download instances into data_dir() and return their paths.

    download('rt-qp-0001')          one instance
    download(['rt-qp-0001', ...])   a list of names, e.g. the result of select(...)
    download('qp')                  every instance of a problem class
    download('qp', 'lung', type=1)  filtered by disease site / any manifest column
    Files already present with the right size and md5 are skipped."""
    names = {r['name'] for r in list_instances()}
    if isinstance(name_or_class, (list, tuple, set)):
        if disease_site is not None or filters:
            raise ValueError('filters only apply when downloading a problem class')
        unknown = sorted(set(name_or_class) - names)
        if unknown:
            raise KeyError('unknown instance(s): {}'.format(', '.join(unknown)))
        wanted = list(dict.fromkeys(name_or_class))             # keep order, drop repeats
    elif name_or_class in names:
        if disease_site is not None or filters:
            raise ValueError('filters only apply when downloading a problem class')
        wanted = [name_or_class]
    elif name_or_class.lower() in PROBLEM_CLASSES:
        if disease_site is not None:
            filters['disease_site'] = disease_site
        wanted = select(problem_class=name_or_class, **filters)
    else:
        raise KeyError('not an instance name or problem class: {}'.format(name_or_class))
    rows = [info(n) for n in wanted]
    todo = [r for r in rows if not _is_cached(data_dir() / (r['name'] + '.h5'), r)]
    print('{} instance(s), {} to download ({:.2f} GB)'.format(
        len(rows), len(todo), sum(r['file_size_bytes'] for r in todo) / 1e9))
    for i, r in enumerate(todo, start=1):
        if not r.get('url'):
            raise RuntimeError('{} has no download url yet'.format(r['name']))
        dst = data_dir() / (r['name'] + '.h5')
        part = dst.with_suffix('.h5.part')
        print('[{}/{}] {} ({:.2f} GB)'.format(i, len(todo), r['name'], r['file_size_bytes'] / 1e9), flush=True)
        with urllib.request.urlopen(r['url'], timeout=60) as resp, open(part, 'wb') as out:
            shutil.copyfileobj(resp, out, length=1 << 24)
        if r.get('md5') and _md5(part) != r['md5']:
            part.unlink()
            raise RuntimeError('md5 mismatch for {} -- download removed, please retry'.format(r['name']))
        os.replace(part, dst)
    return [data_dir() / (r['name'] + '.h5') for r in rows]


class Problem:
    """A QP  minimize 0.5 x'Qx + q'x + c  subject to  Gx <= h.

    Q and G are scipy CSR matrices, q and h numpy vectors, c a float. `info` holds the key/value
    metadata stored in the file (name, disease_site, modality, type, protocol, ...)."""

    def __init__(self, name, Q, q, c, G, h, info, format_note):
        self.name, self.Q, self.q, self.c, self.G, self.h = name, Q, q, c, G, h
        self.info, self.format_note = info, format_note

    @property
    def n(self):
        return self.Q.shape[0]

    @property
    def m(self):
        return self.G.shape[0]

    def objective(self, x):
        x = np.asarray(x, dtype=np.float64)
        return float(0.5 * x @ (self.Q @ x) + self.q @ x + self.c)

    def __repr__(self):
        return '<Problem {}: n={} m={} nnz_Q={} nnz_G={}>'.format(self.name, self.n, self.m, self.Q.nnz, self.G.nnz)


def _read_str(f, key):
    v = f[key][()]
    return v.decode('utf-8') if isinstance(v, bytes) else str(v)


def _read_csr(f, key):
    return sp.csr_matrix((f[key + '_data'][()], f[key + '_indices'][()], f[key + '_indptr'][()]),
                         shape=tuple(int(s) for s in f[key + '_shape'][()]))


def load(name_or_path):
    """Load an instance by name (downloaded on first use) or from a local .h5 path."""
    path = Path(name_or_path)
    if not (path.suffix == '.h5' and path.exists()):
        path = download(str(name_or_path))[0]
    with h5py.File(path, 'r') as f:
        meta = {}
        for line in _read_str(f, 'data_note').splitlines():
            if ':' in line:
                k, v = line.split(':', 1)
                meta[k.strip()] = v.strip()
        return Problem(name=meta.get('name', path.stem), Q=_read_csr(f, 'Q'), q=f['q'][()], c=float(f['c'][()]),
                       G=_read_csr(f, 'G'), h=f['h'][()], info=meta, format_note=_read_str(f, 'format_note'))


def load_solution(name_or_path, solver='warmip'):
    """Load a reference solution: a local <name>-sol-<solver>.h5 path, or an instance name looked up as
    data_dir()/<name>-sol-<solver>.h5 (solver: lowercase name, default 'warmip').

    Returns a dict with name, x (length n), dual (length m, multipliers of Gx <= h), objective,
    primal_res, dual_res, gap, dual_infeas, solver and format_note. Reference solutions come from a
    numerical solver and are NOT guaranteed to be optimal; check one with
    verify_solution(load(name), sol['x'], sol['dual'])."""
    path = Path(name_or_path)
    if not (path.suffix == '.h5' and path.exists()):
        path = data_dir() / '{}-sol-{}.h5'.format(name_or_path, solver)
        if not path.exists():
            raise FileNotFoundError('no {} solution file for {} (looked for {})'.format(solver, name_or_path, path))
    with h5py.File(path, 'r') as f:
        ref = {'name': path.stem.split('-sol-')[0], 'x': f['x'][()], 'dual': f['dual'][()],
               'solver': _read_str(f, 'solver'), 'format_note': _read_str(f, 'format_note')}
        for k in ('objective', 'primal_res', 'dual_res', 'gap', 'dual_infeas'):
            ref[k] = float(f[k][()])
    return ref