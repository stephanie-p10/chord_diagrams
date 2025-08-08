"""This module contains the algorithm that expands non chord patterns into chord patterns"""
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))


from typing import TYPE_CHECKING, Iterable, Tuple, List

from itertools import product, chain

from chords import Chord, GriddedChord
from tiling import Tiling
Cell = Tuple[int, int]

class Expansion:
    """
    This class contains method for expanding the non chord patterns 
    in obstructions and requirements to chord patterns
    """
    def __init__(
            self,
            obstructions: List["GriddedChord"],
            requirements: List[Tuple["GriddedChord", ...]],
            dimensions: Tuple[int, int],
            cells: Tuple[Tuple[int, int]] = None
    ):
        self._obstructions = obstructions
        self._requirements = requirements
        self._dimensions = dimensions

        if cells == None:
            self._cells = list(product(range(0, self._dimensions[0]), range(0, self._dimensions[1])))
        else:
            self._cells = cells
    

    def _expand_single_point(self, gc: GriddedChord, obs) -> Iterable[GriddedChord]:
        """Returns all possible ways a match to the first unpaired chord can be added, with a position in cells.

        If gc has no unpaired points, return none. 
        If it is not possible to expand the first unpaired point of gc into cells, return []"""
        
        gc_pattern = list(gc._patt)
        gc_cells = list(gc._pos)

        chord_to_add = None
        # chord_dict was created as chord: (source idx, sink idx), where sink is -1 if not found
        for chord, indices in sorted(gc.chord_dict.items()): 
            if indices[1] == -1:
                chord_to_add = chord
                break

        if chord_to_add == None: # no unpaired chords were found
            return []
        
        # finds the smallest possible index the point to be added can be inserted into
        smallest_idx = 0
        if chord_to_add > 0:
            smallest_idx = gc.chord_dict[chord_to_add - 1][0] + 1

        # finds the largest possible index the point being added can have
        largest_idx = len(gc_pattern) + 1
        for chord_val in gc.chord_dict.keys():
            # if the source of any chord larger then chord_to_add comes before the point of chord_to_add,
            # the endpoint being placed in chord_to_add must come before this source.
            if chord_val > chord_to_add and gc.chord_dict[chord_val][0] < gc.chord_dict[chord_to_add][0]:
                largest_idx = gc.chord_dict[chord_val][0] + 1

        expanded_chords = []

        # loops over all possible indices where the chord being added could be inserted
        for idx in range(smallest_idx, largest_idx):
            possible_new_patt = gc_pattern.copy()
            possible_new_patt.insert(idx, chord_to_add)

            # loops over all possible cells the new chord could be added into
            for cell in self._cells:
                possible_new_poslist = gc_cells.copy()
                possible_new_poslist.insert(idx, cell)
                new_gc = GriddedChord(Chord(possible_new_patt), possible_new_poslist)
                # if the expansion is consistant with the tiling, add it to the ways that this chord can be expanded
                if not new_gc.contradictory() and new_gc.avoids(*obs):
                    expanded_chords.append(new_gc)

        return list(set(expanded_chords))
        
    def expand_gridded_chord(self, gc: GriddedChord, obs) -> Iterable[GriddedChord]:
        chord = gc._chord
        if chord.is_valid_chord():
            return [gc]
        
        chords_to_build = [gc]

        count = 0
        while (chords_to_build != [] and not chords_to_build[0]._chord.is_valid_chord()):
            count += 1
            extened_chords_to_build = []
            for chord in chords_to_build:
                extened_chords_to_build += self._expand_single_point(chord, obs)
            chords_to_build = list(set(extened_chords_to_build))
        
        return chords_to_build

    def expand_gridded_chords(self, 
                              gcs_to_expand: Iterable["GriddedChord"],
                              gcs_to_avoid: Iterable["GriddedChord"]):
        new_gcs = []
        gcs = sorted(gcs_to_expand)
        for gc in gcs:
            new_gcs += self.expand_gridded_chord(gc, gcs_to_avoid)

        return new_gcs

    def expand_obstructions(self) -> None:
        obs = self._obstructions
        new_obs = []

        for idx, ob in enumerate(obs):
            # we want to delete and then expand this 
            expanded = self.expand_gridded_chord(ob, new_obs)
            for new_ob in expanded:
                new_obs.append(new_ob)

        self._obstructions = new_obs

    # sToDo: this is slightly inefficient, where since we expand requirments regarless if they contain smaller ones in 
    # the same list. If this is the case, the bigger requirement is redundant
    def expand_requirements(self):
        self._requirements = list((self.expand_gridded_chords(reqlist, self._obstructions) for reqlist in self._requirements))

