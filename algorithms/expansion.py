"""This module contains the algorithm that expands non chord patterns into chord patterns"""
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))


from typing import TYPE_CHECKING, Iterable

from chords import Chord, GriddedChord
from tiling import Tiling

class Expansion:
    """
    This class contains method for expanding the non chord patterns 
    in obstructions and requirements to chord patterns
    """
    def __init__(
            self,
            obstructions: tuple["GriddedChord", ...],
            requirements: tuple[tuple["GriddedChord", ...], ...],
            dimensions: tuple[int, int],
    ):
        self._obstructions = obstructions
        self._requirements = requirements
        self._dimensions = dimensions

    def expand_gridded_chord(self, gc: GriddedChord) -> Iterable[GriddedChord]:
        chord = gc._chord
        if chord.is_valid_chord():
            return (gc)
        
        else:
            
            return

    def expand_gridded_chords(self, gcs: Iterable["GriddedChord"]):
        new_gcs = []
        for gc in gcs:
            new_gcs += self.expand_gridded_chord(gc)

        return new_gcs

    def expand_obstructions(self):
        self._obstructions = self.expand_gridded_chords(self._obstructions)

    def expand_requirements(self):
        self._requirements = ((self.expand_gridded_chords(req) for req in reqlist) for reqlist in self._requirements)

    

