"""Make `src/` imports work when running example files directly."""

import sys
from pathlib import Path


repo_root = Path(__file__).resolve().parent.parent
src = repo_root / "src"
if src.is_dir():
    sys.path.insert(0, str(src))

