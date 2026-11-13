"""radoptlib2mps: convert RadOptLIB h5 instances to MPS (QPS) files.

    radoptlib2mps FOLDER             every instance in a folder
    radoptlib2mps rt-qp-0001.h5      single files
    radoptlib2mps                    every instance in the current folder

The MPS files are written next to the h5 files (or into --out) as <name>.mps.gz. The conversion itself is
radoptlib.to_mps; this module only finds the files, names the output and reports the result.
"""
import argparse
import os
import sys
import time
from pathlib import Path

import h5py

REQUIRED = ('Q_data', 'Q_indices', 'Q_indptr', 'Q_shape', 'G_data', 'G_indices', 'G_indptr', 'G_shape', 'q', 'h', 'c')
EXAMPLES = """examples:
  radoptlib2mps data                 convert every instance in the folder "data"
  radoptlib2mps data --out mps       ... and write the MPS files into the folder "mps"
  radoptlib2mps rt-qp-0001.h5        convert one file
  radoptlib2mps                      convert every instance in the current folder
"""


class CliError(Exception):
    """A problem with the command line (exit code 2)."""


def _is_solution(path):
    return '-sol-' in path.name


def _is_instance(path):
    """True if the file is an h5 file with the datasets of an instance."""
    try:
        with h5py.File(path, 'r') as f:
            return all(k in f for k in REQUIRED)
    except OSError:
        return False


def _collect(paths):
    """Expand the command-line paths into (files, notes): the h5 files to convert and one line per file left out."""
    files, notes = [], []
    for p in (Path(x) for x in paths):
        if p.is_dir():
            found = sorted(f for f in p.iterdir() if f.is_file() and f.suffix.lower() == '.h5')
            if not found:
                notes.append('{}: no .h5 files in this folder'.format(p))
            files.extend(f for f in found if not _is_solution(f))        # solution files are left out silently
        elif p.is_file():
            if p.suffix.lower() != '.h5':
                raise CliError('{} is not an .h5 file'.format(p))
            if _is_solution(p):
                raise CliError('{} is a solution file, not an instance'.format(p))
            files.append(p)
        else:
            raise CliError('{} does not exist'.format(p))
    seen, unique = set(), []
    for f in files:
        key = f.resolve()
        if key not in seen:
            seen.add(key)
            unique.append(f)
    return unique, notes


def _size(n_bytes):
    return '{:.2f} GB'.format(n_bytes / 1e9) if n_bytes >= 1e9 else '{:.1f} MB'.format(n_bytes / 1e6)


def _convert_one(src, dst):
    """Convert one file. The MPS file is written under a temporary name and renamed when it is complete,
    so an interrupted run never leaves a partial file under the final name."""
    from . import load, to_mps                      # imported here so that --help and --version stay fast
    tmp = dst.with_name('partial-' + dst.name)       # keeps the .mps/.mps.gz ending, which selects the compression
    try:
        problem = load(src)
        to_mps(problem, tmp)
        os.replace(tmp, dst)
    finally:
        if tmp.exists():
            tmp.unlink()
    return problem


def _build_parser():
    from . import __version__
    p = argparse.ArgumentParser(
        prog='radoptlib2mps', formatter_class=argparse.RawDescriptionHelpFormatter, epilog=EXAMPLES,
        description='Convert RadOptLIB instances (.h5) to MPS files (free-format QPS: MPS with a QUADOBJ section).')
    p.add_argument('paths', nargs='*', metavar='FOLDER_OR_FILE',
                   help='folders with .h5 files and/or single .h5 files (default: the current folder)')
    p.add_argument('--out', metavar='FOLDER', help='write the MPS files into this folder (default: next to each .h5 file)')
    p.add_argument('--no-gzip', action='store_true', help='write plain .mps files instead of .mps.gz (5 to 6 times larger)')
    p.add_argument('--force', action='store_true', help='convert again and overwrite MPS files that already exist')
    p.add_argument('--version', action='version', version='radoptlib2mps ' + __version__)
    return p


def run(argv=None):
    args = _build_parser().parse_args(argv)
    files, notes = _collect(args.paths or ['.'])
    for line in notes:
        print(line)
    if not files:
        where = 'the current folder' if not args.paths else 'the given location'
        raise CliError('no instance files (.h5) found in {}. Run "radoptlib2mps --help" for the usage.'.format(where))
    out_dir = Path(args.out) if args.out else None
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
    ending = '.mps' if args.no_gzip else '.mps.gz'
    converted, existing, skipped, failed = [], [], [], []
    for i, src in enumerate(files, start=1):
        dst = (out_dir or src.parent) / (src.stem + ending)
        tag = '[{}/{}] {}'.format(i, len(files), src.name)
        if not _is_instance(src):
            print('{}: skipped, not a RadOptLIB instance'.format(tag))
            skipped.append(src)
            continue
        if dst.exists() and not args.force:
            print('{}: {} already exists (use --force to convert again)'.format(tag, dst.name))
            existing.append(src)
            continue
        print('{}: converting ({}) ...'.format(tag, _size(src.stat().st_size)), flush=True)
        start = time.time()
        try:
            problem = _convert_one(src, dst)
        except KeyboardInterrupt:
            print('interrupted; {} was not written'.format(dst.name))
            raise
        except MemoryError:
            print('    FAILED: not enough memory for this instance')
            failed.append(src)
            continue
        except Exception as e:                       # keep going with the remaining files
            print('    FAILED: {}: {}'.format(type(e).__name__, e))
            failed.append(src)
            continue
        print('    {} variables, {} constraints -> {} ({}, {:.0f} s)'.format(
            problem.n, problem.m, dst.name, _size(dst.stat().st_size), time.time() - start), flush=True)
        converted.append(src)
    parts = ['{} converted'.format(len(converted))]
    if existing:
        parts.append('{} already there'.format(len(existing)))
    if skipped:
        parts.append('{} skipped'.format(len(skipped)))
    if failed:
        parts.append('{} failed ({})'.format(len(failed), ', '.join(f.name for f in failed)))
    print('Done: ' + ', '.join(parts) + '.')
    return 1 if failed else 0


def _pause_if_double_clicked():
    """Started by double-click on Windows, the console window closes as soon as the program ends.
    If this program owns the console (no shell is attached to it), wait for Enter so the result can be read."""
    if os.name != 'nt' or not getattr(sys, 'frozen', False):
        return
    try:
        import ctypes
        buf = (ctypes.c_uint * 4)()
        attached = ctypes.windll.kernel32.GetConsoleProcessList(buf, 4)
        if 0 < attached <= 2:                        # the program and its unpacking launcher only
            input('Press Enter to close this window.')
    except Exception:
        pass


def main(argv=None):
    try:
        code = run(argv)
    except CliError as e:
        print('radoptlib2mps: ' + str(e), file=sys.stderr)
        code = 2
    except KeyboardInterrupt:
        code = 130
    _pause_if_double_clicked()
    return code


if __name__ == '__main__':
    sys.exit(main())
