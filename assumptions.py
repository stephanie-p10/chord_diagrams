import abc
from importlib import import_module
from itertools import chain
from typing import TYPE_CHECKING, Iterable, List, Optional, Tuple, Type

from permuta import Perm

from chords import GriddedChord

Cell = Tuple[int, int]

if TYPE_CHECKING:
    from .tiling import Tiling


class TrackingAssumption:
    """
    An assumption used to keep track of the occurrences of a set of gridded
    chord diagrams.
    """

    def __init__(self, gcs: Iterable[GriddedChord]):
        self.gcs = tuple(sorted(set(gcs)))

    @classmethod
    def from_cells(cls, cells: Iterable[Cell]) -> "TrackingAssumption":
        gcs = [GriddedChord.single_cell((0,0), cell) for cell in cells]
        return cls(gcs)

    def get_cells(self) -> Tuple[Cell, ...]:
        assert all(len(gc) == 1 for gc in self.gcs)
        return tuple(gc.pos[0] for gc in self.gcs)

    def avoiding(
        self,
        obstructions: Iterable[GriddedChord],
        active_cells: Optional[Iterable[Cell]] = None,
    ) -> "TrackingAssumption":
        """
        Return the tracking absumption where all of the gridded chords avoiding
        the obstructions are removed. If active_cells is not None, then only assumputions 
        completely contained in active cells will remain.
        """
        obstructions = tuple(obstructions)
        if active_cells is not None:
            return self.__class__(
                tuple(
                    gc
                    for gc in self.gcs
                    if all(cell in active_cells for cell in gc.pos)
                    and gc.avoids(*obstructions)
                )
            )
        return self.__class__(tuple(gc for gc in self.gcs if gc.avoids(*obstructions)))

    def get_value(self, gc: GriddedChord) -> int:
        """
        Return the number of occurrences of each of the gridded chord diagrams being tracked in
        the gridded chord diagram gc.
        """
        return len(list(chain.from_iterable(p.occurrences_in(gc) for p in self.gcs)))

#    def get_components(self, tiling: "Tiling") -> List[List[GriddedChord]]:
#        """
#        Return the lists of gcs that count exactly one occurrence.
#        Only implemented for when a size one gc is in a point cell.
#        """
#        return [
#            [gp] for gp in self.gcs if len(gp) == 1 and gp.pos[0] in tiling.point_cells
#        ]
#
#    def remove_components(self, tiling: "Tiling") -> "TrackingAssumption":
#        """
#        Return the TrackingAssumption found by removing all the components
#        found by the get_components method.
#        """
#        gps_to_remove = set(chain.from_iterable(self.get_components(tiling)))
#        return self.__class__(gp for gp in self.gcs if gp not in gps_to_remove)

    def to_jsonable(self) -> dict:
        """Return a dictionary form of the assumption."""
        c = self.__class__
        return {
            "class_module": c.__module__,
            "assumption": c.__name__,
            "gcs": [gc.to_jsonable() for gc in self.gcs],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "TrackingAssumption":
        """Return the assumption from the json dict representation."""
        module = import_module(d["class_module"])
        AssumpClass: Type["TrackingAssumption"] = getattr(module, d["assumption"])
        assert issubclass(
            AssumpClass, TrackingAssumption
        ), "Not a valid TrackingAssumption"
        gcs = [GriddedChord.from_dict(gc) for gc in d["gcs"]]
        return AssumpClass(gcs)

    def __eq__(self, other) -> bool:
        if other.__class__ == TrackingAssumption:
            return bool(self.gcs == other.gcs)
        return NotImplemented

    def __lt__(self, other) -> bool:
        if isinstance(other, TrackingAssumption):
            key_self = (self.__class__.__name__, self.gcs)
            key_other = (other.__class__.__name__, other.gcs)
            return key_self < key_other
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.gcs)

    def __repr__(self) -> str:
        return self.__class__.__name__ + f"({self.gcs})"

    def __str__(self):
        if all(len(gp) == 1 for gp in self.gcs):
            cells = ", ".join(str(gc.pos[0]) for gc in self.gcs)
            return f"can count points in cell{'s' if len(self.gcs) > 1 else ''} {cells}"
        return f"can count occurrences of {', '.join(str(gc) for gc in self.gcs)}"