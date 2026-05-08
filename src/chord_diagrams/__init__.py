"""Top-level package for the chord-diagram tilings library.

Most code should import the core primitives from here:

- `Chord`: chord-diagram pattern
- `GriddedChord`: a chord pattern embedded into grid cells
"""

from .chords import Chord, GriddedChord

__version__ = "4.0.0"

__all__ = ["Chord", "GriddedChord"]
