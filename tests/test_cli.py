"""radoptlib2mps command line: file discovery, output names, skipping and error handling."""
import gzip

import h5py
import pytest

import radoptlib
from radoptlib import cli

from conftest import ATTRS, make_qp, write_h5


@pytest.fixture
def folder(tmp_path):
    """Two instances, one reference solution file and one unrelated h5 file."""
    d = tmp_path / 'data'
    d.mkdir()
    write_h5(d / 'rt-qp-0001.h5', make_qp(seed=1), attrs=dict(ATTRS, name='rt-qp-0001'))
    write_h5(d / 'rt-qp-0002.h5', make_qp(seed=2), attrs=dict(ATTRS, name='rt-qp-0002'))
    with h5py.File(d / 'rt-qp-0001-sol-warmip.h5', 'w') as f:
        f.create_dataset('x', data=[0.0])
    with h5py.File(d / 'other.h5', 'w') as f:
        f.create_dataset('something', data=[1.0])
    return d


def test_folder(folder, capsys):
    assert cli.main([str(folder)]) == 0
    assert sorted(p.name for p in folder.glob('*.mps*')) == ['rt-qp-0001.mps.gz', 'rt-qp-0002.mps.gz']
    out = capsys.readouterr().out
    assert 'other.h5: skipped, not a RadOptLIB instance' in out and 'sol-warmip' not in out
    assert 'Done: 2 converted, 1 skipped.' in out
    # same content as the library function
    radoptlib.to_mps(radoptlib.load(folder / 'rt-qp-0001.h5'), folder / 'ref.mps.gz')
    assert gzip.open(folder / 'rt-qp-0001.mps.gz').read() == gzip.open(folder / 'ref.mps.gz').read()


def test_current_folder_is_the_default(folder, monkeypatch):
    monkeypatch.chdir(folder)
    assert cli.main([]) == 0
    assert (folder / 'rt-qp-0001.mps.gz').exists() and (folder / 'rt-qp-0002.mps.gz').exists()


def test_existing_files_are_kept_unless_forced(folder, capsys):
    target = folder / 'rt-qp-0001.mps.gz'
    target.write_bytes(b'old')
    assert cli.main([str(folder)]) == 0
    assert target.read_bytes() == b'old'
    assert '1 converted, 1 already there' in capsys.readouterr().out
    assert cli.main([str(folder), '--force']) == 0
    assert gzip.open(target).read().startswith(b'NAME rt-qp-0001')


def test_single_files_out_folder_and_plain(folder, tmp_path):
    out = tmp_path / 'mps'
    assert cli.main([str(folder / 'rt-qp-0002.h5'), '--out', str(out), '--no-gzip']) == 0
    assert [p.name for p in out.iterdir()] == ['rt-qp-0002.mps']
    assert (out / 'rt-qp-0002.mps').read_text().startswith('NAME rt-qp-0002')


def test_failure_continues_and_leaves_no_partial_file(folder, monkeypatch, capsys):
    real = radoptlib.to_mps

    def flaky(problem, path=None, compress=True):
        real(problem, path, compress)
        if problem.name == 'rt-qp-0001':
            raise RuntimeError('disk full')

    monkeypatch.setattr(radoptlib, 'to_mps', flaky)
    assert cli.main([str(folder)]) == 1
    assert sorted(p.name for p in folder.glob('*.mps*')) == ['rt-qp-0002.mps.gz']
    out = capsys.readouterr().out
    assert 'FAILED: RuntimeError: disk full' in out and '1 converted, 1 skipped, 1 failed (rt-qp-0001.h5)' in out


def test_usage_errors(folder, tmp_path, capsys):
    empty = tmp_path / 'empty'
    empty.mkdir()
    assert cli.main([str(empty)]) == 2
    assert cli.main([str(tmp_path / 'missing')]) == 2
    assert cli.main([str(folder / 'rt-qp-0001-sol-warmip.h5')]) == 2
    err = capsys.readouterr().err
    assert 'no instance files' in err and 'does not exist' in err and 'is a solution file' in err
