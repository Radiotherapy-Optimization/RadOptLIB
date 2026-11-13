# MPS files

RadOptLIB instances are distributed as HDF5 files (`.h5`). For solvers that read MPS files, the
program `radoptlib2mps` converts the downloaded `.h5` files on your own computer. It is a single
file and needs no Python and no installation.

## 1. Get the program

Download the file for your system from the
[latest release](https://github.com/Radiotherapy-Optimization/RadOptLIB/releases/latest):

| System | File |
|---|---|
| Windows (64-bit) | `radoptlib2mps-windows-x64.exe` |
| Linux (64-bit) | `radoptlib2mps-linux-x64` |
| macOS (Apple silicon) | `radoptlib2mps-macos-arm64` |

The programs are not signed, so the system asks for confirmation the first time:

- **Windows:** if "Windows protected your PC" appears, choose "More info" and then "Run anyway".
- **macOS:** in a terminal, run `chmod +x radoptlib2mps-macos-arm64` and
  `xattr -d com.apple.quarantine radoptlib2mps-macos-arm64` once.
- **Linux:** run `chmod +x radoptlib2mps-linux-x64` once.

## 2. Download instances

Download the `.h5` files you want into one folder, for example with the
[dashboard](https://radiotherapy-optimization.github.io/RadOptLIB/dashboard/) and the download script it provides.

## 3. Convert

Give the program the folder that holds the `.h5` files:

```
radoptlib2mps-windows-x64.exe C:\data\radoptlib        (Windows)
./radoptlib2mps-linux-x64 ~/data/radoptlib             (Linux)
./radoptlib2mps-macos-arm64 ~/data/radoptlib           (macOS)
```

Every instance in the folder is converted, and `rt-qp-0001.mps.gz` is written next to
`rt-qp-0001.h5`. Without a folder, the program converts the folder it is started in, so on Windows
it can also be placed next to the `.h5` files and started with a double-click.

| Usage | Effect |
|---|---|
| `radoptlib2mps FOLDER` | convert every instance in the folder (subfolders are not searched) |
| `radoptlib2mps FILE.h5 ...` | convert single files |
| `--out FOLDER` | write the MPS files into another folder |
| `--no-gzip` | write plain `.mps` files instead of `.mps.gz` |
| `--force` | convert again and overwrite MPS files that already exist |
| `--help`, `--version` | show the usage or the version |

MPS files that already exist are kept, so an interrupted run can simply be started again.
Reference solution files (`...-sol-<solver>.h5`) and other `.h5` files in the folder are left alone.
If an instance cannot be converted, the program reports it, continues with the next one and lists
the failures at the end.

## The MPS files

The files are free-format MPS files with a `QUADOBJ` section (also called QPS), the form most QP
solvers read:

- rows: the objective `obj` and one `L` row per row of G, named `c1` ... `cm` in the order of the
  HDF5 file;
- columns: the variables `x1` ... `xn`;
- `BOUNDS`: the bounds l and u of every variable, written explicitly (free variables are marked
  as free);
- `RHS`: the vector h, and the objective constant c stored as `-c` on the objective row;
- `QUADOBJ`: the lower triangle of Q (readers apply the factor ½);
- numbers are written with full double precision, so the MPS file describes exactly the same
  problem as the HDF5 file.

Most solvers read `.mps.gz` files directly. If yours does not, use `--no-gzip` or unpack the file.

## Disk space, time and memory

- A compressed MPS file is about 1.4 to 2 times the size of its `.h5` file; a plain one is 5 to 6
  times the size.
- Time and memory grow with the number of nonzeros of the instance. An instance with 15 million
  nonzeros takes about one minute and 1 GB of memory. The largest instances (more than 400 million
  nonzeros) are estimated to need 15 to 20 GB of memory and half an hour or more each.

## With Python

Installing the `radoptlib` package provides the same command, and the conversion is also available
as a function:

```python
import radoptlib
radoptlib.to_mps(radoptlib.load('rt-qp-0001.h5'), 'rt-qp-0001.mps.gz')
```

## Building the program

The programs on the release page are built from this repository by the workflow
`.github/workflows/build-executables.yml` every time a release is published, so they always contain
the converter of that release. To build one by hand:

```
pip install . pyinstaller
pyinstaller mps/radoptlib2mps.spec
```

The result is `dist/radoptlib2mps` (`dist\radoptlib2mps.exe` on Windows).
