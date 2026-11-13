# PyInstaller recipe for the stand-alone radoptlib2mps program (one file, no Python needed to run it).
# Build from the repository root:   pyinstaller mps/radoptlib2mps.spec
from pathlib import Path

root = Path(SPECPATH).parent

a = Analysis(
    [str(root / 'mps' / 'radoptlib2mps.py')],
    pathex=[str(root / 'src')],
    excludes=['tkinter', 'matplotlib', 'pytest', 'IPython'],
)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, a.binaries, a.datas, [], name='radoptlib2mps', console=True, upx=False)
