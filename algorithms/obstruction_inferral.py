import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

import abc
from typing import TYPE_CHECKING, Iterable, List, Optional, Set, Tuple

from chords import GriddedChord
#from tilings.algorithms.gridded_perm_generation import GriddedPermsOnTiling

if TYPE_CHECKING:
    from tiling import Tiling

Cell = Tuple[int, int]


class ObstructionInferral(abc.ABC):
    """
    Algorithm to compute the tiling created by adding all obstructions from
    `self.potential_new_obs()` that can be added to the tiling.
    """

    def __init__(self, tiling: "Tiling"):
        self._tiling = tiling
        self._new_obs: Optional[List[GriddedChord]] = None

    @abc.abstractmethod
    def potential_new_obs(self) -> Iterable[GriddedChord]:
        """
        Return an iterable of new obstructions that should be added to the
        tiling if possible.
        """
    def new_obs(self, yield_non_minimal: bool = False) -> List[GriddedChord]:
        """
        Returns the list of new obstructions that can be added to the tiling.
        """
        if self._new_obs is not None:
            return self._new_obs

        chords_to_check = tuple(self.potential_new_obs())
        if not chords_to_check: # if there are no chords to check, there are no possible obstructions can add with the given parameters
            self._new_obs = []
            return self._new_obs

        # finds max number of chords a single chord diagram has in the potential new obstruction patterns
        max_len_in_chords_to_check = max(map(lambda chord: max(chord._patt) + 1 if chord else 0, chords_to_check)) 
        #print(max_len_in_chords_to_check)
        # Calculates max length a smallest chord diagram containing a possible new obstruction could have
        max_length = (
            self._tiling.maximum_length_of_minimum_gridded_chord()
            + max_len_in_chords_to_check
        )
        #print(max_len_in_chords_to_check, self._tiling.maximum_length_of_minimum_gridded_chord())
        # A list of all gridded chord diagrams griddable on the tiling of up to max_length 
        # (*2 since the method works in number of points not number of chords)
        chords_on_tiling = self._tiling.all_chords_on_tiling(2 * max_length)

        chords_left = set(chords_to_check) # found by .potential_new_obs()
        #print("chords left generated")
        for gc in chords_on_tiling: # loop over every gridded chord that can be gridded on the tiling, up to max_length
            to_remove: List[GriddedChord] = []
            for chord in chords_left: # for every potential new ob to check
                if gc.contains(chord): # if a possible gridded chord contains this ob
                    to_remove.append(chord) # we need to remove it from the possible obstructions to add
            chords_left.difference_update(to_remove) # updates perms_left so that none conflict with gp
            if not chords_left: # if there are no perms_left we want to stop
                break

        
        self._new_obs = sorted(chords_left) # the remaining chords after they all have been checked.
        return self._new_obs 

    @staticmethod
    def can_add_obstruction(obstruction: GriddedChord, tiling: "Tiling") -> bool:
        """Return true if `obstruction` can be added to `tiling`."""
        return tiling.add_requirement(obstruction.patt, obstruction.pos).is_empty()

    def obstruction_inferral(self) -> "Tiling":
        """
        Return the tiling with the new obstructions.
        """
        return self._tiling.add_obstructions(self.new_obs())

    # TODO: move to strategy class

    def formal_step(self):
        """Return a string describing the operation performed."""
        return f"Added the obstructions {self.new_obs()}."


class SubobstructionInferral(ObstructionInferral):
    """
    Algorithm to compute the tiling created by adding all
    subobstructions which can be added.
    """

    def potential_new_obs(self) -> Set[GriddedChord]:
        """
        Return the set of all subobstructions of the tiling.
        """
        subobs: Set[GriddedChord] = set()
        for ob in self._tiling.obstructions:
            subobs.update(ob.all_subchords(proper=True, return_all_subpatts=True))
        subobs.remove(GriddedChord.empty_chord())
        return subobs


class AllObstructionInferral(ObstructionInferral):
    """
    Algorithm to compute the tiling created by adding all
    obstructions with up to obstruction_length number of points.
    """

    def __init__(self, tiling: "Tiling", obstruction_length: Optional[int]) -> None:
        super().__init__(tiling)
        self._obs_len = obstruction_length

    @property
    def obstruction_length(self) -> Optional[int]:
        return self._obs_len

    def not_required(self, gc: GriddedChord) -> bool:
        """
        Returns True if the gridded chord `gc` is not required by any
        requirement list of the tiling.
        """
        return all(
            any(gc not in req for req in req_list)
            for req_list in self._tiling.requirements
        )

    def potential_new_obs(self) -> List[GriddedChord]:
        """
        Iterator over all possible obstruction of `self.obstruction_length`.
        """
        if not self._tiling.requirements:
            return []
        no_req_tiling = self._tiling.__class__(self._tiling.obstructions)
        n = self._obs_len
        pot_obs = filter(self.not_required, no_req_tiling.all_chords_on_tiling(n, True)) 
        return list(gc for gc in pot_obs)


class EmptyCellInferral(AllObstructionInferral):
    """
    Try to add a point obstruction to all the active non positive cell
    """

    def __init__(self, tiling: "Tiling"):
        super().__init__(tiling, 1)

    def empty_cells(self) -> List[Cell]:
        """
        Return an iterator over all cell that where discovered to be empty.
        """
        return list(ob.pos[0] for ob in self.new_obs())
