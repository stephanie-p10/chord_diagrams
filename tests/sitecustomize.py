"""Make `src/` imports work when running test files directly.

When you run `python tests/test_foo.py`, Python puts `tests/` on `sys.path` but
*not* the repository root or `src/`. Importing this module (which Python does
automatically if it is on `sys.path`) fixes that by adding `src/` to `sys.path`.
"""

import sys
from pathlib import Path


def _add_src_to_path() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    src = repo_root / "src"
    if src.is_dir():
        sys.path.insert(0, str(src))


_add_src_to_path()

