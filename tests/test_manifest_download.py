"""Manifest lookup and downloads against a fake manifest whose urls are local file:// paths."""
import hashlib

import pytest

import radoptlib
from radoptlib import io as rio


@pytest.fixture
def fake_library(tmp_path, monkeypatch, qp_file):
    store = tmp_path / 'store'
    store.mkdir()
    rows = []
    for i, (site, typ) in enumerate([('lung', 0), ('lung', 1), ('prostate', 1)], start=1):
        name = 'rt-qp-{:04d}'.format(i)
        src = store / (name + '.h5')
        src.write_bytes(qp_file.read_bytes())
        rows.append({'name': name, 'problem_class': 'qp', 'problem_type': 'convex continuous QP',
                     'disease_site': site, 'modality': 'vmat', 'type': typ, 'patient': site + '-0001',
                     'protocol': '', 'n_vars': 12, 'n_ineq_cons': 20, 'n_eq_cons': 0, 'nnz_Q': 0, 'nnz_G': 0,
                     'file_size_bytes': src.stat().st_size, 'url': src.as_uri(), 'added_in_version': '1.0',
                     'note': '', 'md5': hashlib.md5(src.read_bytes()).hexdigest()})
    monkeypatch.setattr(rio, '_load_manifest', lambda cls: rows)
    monkeypatch.setenv('RADOPTLIB_DATA', str(tmp_path / 'cache'))
    return rows


def test_select_and_info(fake_library):
    assert radoptlib.select(problem_class='qp') == ['rt-qp-0001', 'rt-qp-0002', 'rt-qp-0003']
    assert radoptlib.select(disease_site='lung', type=1) == ['rt-qp-0002']
    assert radoptlib.select(disease_site=['lung', 'prostate'], type='1') == ['rt-qp-0002', 'rt-qp-0003']
    assert radoptlib.info('rt-qp-0003')['disease_site'] == 'prostate'
    with pytest.raises(KeyError):
        radoptlib.info('rt-qp-9999')
    with pytest.raises(KeyError):
        radoptlib.select(no_such_column=1)


def test_download_variants(fake_library, capsys):
    assert [p.name for p in radoptlib.download('rt-qp-0002')] == ['rt-qp-0002.h5']
    assert [p.name for p in radoptlib.download('qp', 'lung')] == ['rt-qp-0001.h5', 'rt-qp-0002.h5']
    assert len(radoptlib.download('qp')) == 3
    out = capsys.readouterr().out
    assert '1 instance(s), 1 to download' in out
    assert '2 instance(s), 1 to download' in out          # rt-qp-0002 already cached
    assert '3 instance(s), 1 to download' in out
    assert radoptlib.load('rt-qp-0003').n == 12            # load by name uses the cache


def test_select_size_and_ranges(fake_library):
    for r, (n, m, size) in zip(fake_library, [(10, 20, 'small'), (300, 400, 'medium'), (5000, 9000, 'large')]):
        r.update({'n_vars': n, 'n_ineq_cons': m, 'size': size})
    assert radoptlib.select(size='large') == ['rt-qp-0003']
    assert radoptlib.select(size=['small', 'medium']) == ['rt-qp-0001', 'rt-qp-0002']
    assert radoptlib.select(min_n_vars=300) == ['rt-qp-0002', 'rt-qp-0003']
    assert radoptlib.select(min_n_vars=11, max_n_ineq_cons=400) == ['rt-qp-0002']
    assert radoptlib.select(n_eq_cons=0) == ['rt-qp-0001', 'rt-qp-0002', 'rt-qp-0003']
    assert radoptlib.select(max_file_size_bytes=1e12, type=1) == ['rt-qp-0002', 'rt-qp-0003']
    with pytest.raises(KeyError):
        radoptlib.select(max_disease_site=1)                   # not a numeric column
    with pytest.raises(KeyError):
        radoptlib.select(max_nothing=1)


def test_download_selected_list(fake_library):
    selected = radoptlib.select(disease_site=['lung', 'prostate'], type=1)
    assert [p.name for p in radoptlib.download(selected)] == ['rt-qp-0002.h5', 'rt-qp-0003.h5']
    with pytest.raises(KeyError):
        radoptlib.download(['rt-qp-0002', 'rt-qp-9999'])
    with pytest.raises(ValueError):
        radoptlib.download(selected, 'lung')


def test_download_md5_mismatch(fake_library):
    fake_library[0]['md5'] = '0' * 32
    with pytest.raises(RuntimeError, match='md5 mismatch'):
        radoptlib.download('rt-qp-0001')
    assert not (radoptlib.data_dir() / 'rt-qp-0001.h5').exists()


def test_load_solution_downloads_from_manifest(fake_library, qp, tmp_path):
    import h5py
    src = tmp_path / 'store' / 'rt-qp-0001-sol-warmip.h5'
    with h5py.File(src, 'w') as f:
        for k, v in (('x', qp['x']), ('dual_ineq_cons', qp['y']), ('dual_lower_bound', qp['zl']),
                     ('dual_upper_bound', qp['zu'])):
            f.create_dataset(k, data=v)
        f.create_dataset('solver', data='test solver')
        f.create_dataset('format_note', data='reference solution')
    with pytest.raises(FileNotFoundError):                    # no sol_url in the manifest yet
        radoptlib.load_solution('rt-qp-0001')
    fake_library[0].update({'sol_url': src.as_uri(), 'sol_md5': '0' * 32})
    with pytest.raises(RuntimeError, match='md5 mismatch'):
        radoptlib.load_solution('rt-qp-0001')
    assert not (radoptlib.data_dir() / 'rt-qp-0001-sol-warmip.h5').exists()
    fake_library[0]['sol_md5'] = hashlib.md5(src.read_bytes()).hexdigest()
    sol = radoptlib.load_solution('rt-qp-0001')
    assert (sol['x'] == qp['x']).all() and (radoptlib.data_dir() / 'rt-qp-0001-sol-warmip.h5').exists()


def test_download_bad_arguments(fake_library):
    with pytest.raises(KeyError):
        radoptlib.download('lp')
    with pytest.raises(ValueError):
        radoptlib.download('rt-qp-0001', 'lung')