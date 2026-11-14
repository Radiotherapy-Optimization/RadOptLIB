"""Entry point of the stand-alone radoptlib2mps program (built with PyInstaller from radoptlib2mps.spec, see docs/mps.md)."""
import sys

from radoptlib.cli import main

if __name__ == '__main__':
    sys.exit(main())
