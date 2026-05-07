"""Bootstrap so `python3 tests/test_*.py` works with a `src/` layout.

Running a test module directly puts `tests/` on `sys.path` but not `src/`, so
`import chord_diagrams` would fail unless PYTHONPATH is set. Import this module
as the first statement in each test file to make direct execution work.
"""

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir():
    sys.path.insert(0, str(_SRC))

