"""Lightweight compatibility layer for the external `tilings` package.

This repository historically imported helpers from `tilings.*`. In environments where the
external dependency isn't available, we provide the small subset needed by this codebase.
"""

from .exception import InvalidOperationError

__all__ = ["InvalidOperationError"]

