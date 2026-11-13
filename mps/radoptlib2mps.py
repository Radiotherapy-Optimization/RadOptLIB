"""Entry point of the stand-alone radoptlib2mps program (built with PyInstaller, see README.md in this folder)."""
import sys

from radoptlib.cli import main

if __name__ == '__main__':
    sys.exit(main())
