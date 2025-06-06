import json
from itertools import chain, filterfalse, product
from typing import (Any, Callable, Dict, FrozenSet, Iterable, Iterator, List, Optional, Set, Tuple)

from typing_extensions import TypedDict

from comb_spec_searcher import CombinatorialClass, VerificationStrategy
from comb_spec_searcher.exception import StrategyDoesNotApply
from comb_spec_searcher.typing import Parameters

from chords import GriddedChord, Chord
from tilings.misc import intersection_reduce, union_reduce
from tilings.assumptions import (
    ComponentAssumption,
    SkewComponentAssumption,
    SumComponentAssumption,
    TrackingAssumption,
)

__all__ = ["Tiling"]

Cell = Tuple[int, int]
ReqList = Tuple[GriddedChord, ...]
CellBasis = Dict[Cell, Tuple[List[Chord], List[Chord]]]
CellFrozenSet = FrozenSet[Cell]
Dimension = Tuple[int, int]
GCTuple = Tuple[GriddedChord, ...]


# sToDo: make way to handle single point obstructions/requirements

class Tiling(CombinatorialClass):
    """Tiling class.

    Zero-indexed coordinates/cells from bottom left corner where the (x, y)
    cell is the cell in the x-th column and y-th row.

    Tilings store the obstructions, requirements and linkages but also caches the empty
    cells and the active cells.
    """
    def __init__(self,
        dimensions: Iterable[int] = (1, 1),
        obstructions: Iterable[GriddedChord] = tuple(),
        requirements: Iterable[Iterable[GriddedChord]] = tuple(), # might need a case for requirement of having one point in a cell...
        linkages: Iterable[Iterable[Cell]] = tuple(),
        assumptions: Iterable[TrackingAssumption] = tuple()):

        super().__init__()
        self._length = dimensions[0]
        self._height = dimensions[1]
        self._linkages = linkages
        self._obstructions = obstructions
        self._requirements = requirements
        self._assumptions = assumptions

    @property
    def obstructions(self):
        return self._obstructions
    
    @property
    def requirements(self):
        return self._requirements
    
    @property
    def linkages(self):
        return self._linkages

    @property
    def assumptions(self):
        return self._assumptions
    
    def is_empty(self) -> bool:
        """Checks if the tiling is empty.

        Tiling is empty if it has been inferred to be contradictory due to
        contradicting requirements and obstructions or no gridded chord
        can be gridded on the tiling.
        """
        # if any obstruction is empty, no chords can be gridded since all grids contain the empty grid.
        if any(ob.is_empty() for ob in self.obstructions):
            return True
        
        # if any req list is empty or has all empty grids, only the empty grid can be gridded. 
        if any(all(req.is_empty() for req in req_list) for req_list in self.requirements) and len(self.requirements) >= 1:
            return True
        
        sum_max_reqs = 0
        for req_list in self.requirements:
            #if len(req_list) == 0:
                #return True
            req_lengths = [len(req) for req in req_list]
            sum_max_reqs += max(req_lengths)

        if sum_max_reqs == 0: # then there are no requirements, but we still need a chord of size one
            sum_max_reqs += 1
        
        # proved maximum size of smallest chord that can be gridded:
        max_len = sum_max_reqs * 2 - 1
        print(max_len)

        all_chords = []
        for num_chords in range(1, max_len + 1):
            all_chords += list(Chord.of_length(num_chords))

        # product of length and width to get valid cells can probaby be much improved.
        cells = list(product(range(self._length), range(self._height)))
        all_gridded_chords = []
        for chord in all_chords:
            all_gridded_chords += list(GriddedChord.all_grids(chord, cells))

        # now check if any of the chords work??
        #for gc in all_gridded_chords:
        #    if self.contains(gc):
        #        print(gc)
        return not any(self.contains(gc) for gc in all_gridded_chords)

    def contains(self, gc: GriddedChord) -> bool:
        has_reqs = all(gc.contains(*req) for req in self._requirements)
        avoids_ob = not any(gc.contains(ob) for ob in self._obstructions)
        links_connected = all(gc.is_connected(cells) for cells in self._linkages)
        return has_reqs and avoids_ob and links_connected

    def __hash__(self) -> int:
        return (
            hash(self._requirements)
            ^ hash(self._obstructions)
            ^ hash(self._assumptions)
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Tiling):
            return False
        return (
            self.obstructions == other.obstructions
            and self.requirements == other.requirements
            and self.linkages == other.linkages
            and self.assumptions == other.assumptions
        )

    def __ne__(self, other: object) -> bool:
        if not isinstance(other, Tiling):
            return True
        return (
            self.obstructions != other.obstructions
            or self.requirements != other.requirements
            or self.linkages != other.linkages
            or self.assumptions != other.assumptions
        )

    # sToDo: fix for linkages
    def __contains__(self, gp: GriddedChord) -> bool:
        """Test if a gridded chord is griddable on the given tiling."""
        return gp.avoids(*self.obstructions) and all(
            gp.contains(*req) for req in self.requirements
        )

    def __repr__(self) -> str:
        format_string = "Tiling(obstructions={}, requirements={}, assumptions={})"
        #huh I'm not actually sure why this is here. Why don't I want single chords in my obstructions?
        non_point_obstructions = tuple(
            filterfalse(GriddedChord.is_single_chord, self.obstructions)
        )
        return format_string.format(
            non_point_obstructions, self.requirements, self.assumptions
        )

    # not sure this needs to be this long... I don't think I need it printed nice?
    def __str__(self) -> str:
        '''# pylint: disable=too-many-locals
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-statements
        dim_i, dim_j = self.dimensions
        result = []
        # Create tiling lines
        for j in range(2 * dim_j + 1):
            for i in range(2 * dim_i + 1):
                # Whether or not a vertical line and a horizontal line is
                # present
                vertical = i % 2 == 0
                horizontal = j % 2 == 0
                if vertical:
                    if horizontal:
                        result.append("+")
                    else:
                        result.append("|")
                elif horizontal:
                    result.append("-")
                else:
                    result.append(" ")
            result.append("\n")

        labels: Dict[Tuple[Tuple[Chord, ...], bool], str] = {}

        # Put the sets in the tiles

        # How many characters are in a row in the grid
        row_width = 2 * dim_i + 2
        curr_label = 1
        for cell, gridded_chords in sorted(self.cell_basis().items()):
            obstructions, _ = gridded_chords
            basis = list(sorted(obstructions))
            if basis == [Chord((0,0, 1, 1))]:
                continue
            # the block, is the basis and whether or not positive
            block = (tuple(basis), cell in self.positive_cells)
            label = labels.get(block)
            if label is None:
                if basis == [Chord((0, 0, 1, 1)), Chord((0, 1, 0,1))]:
                    if cell in self.positive_cells:
                        label = "\u25cf"
                    else:
                        label = "\u25cb"
                elif basis == [Chord((0, 0, 1, 1))]:
                    label = "\\"
                elif basis == [Chord((0,1,0,1))]:
                    label = "/"
                else:
                    label = str(curr_label)
                    curr_label += 1
                labels[block] = label
            row_index_from_top = dim_j - cell[1] - 1
            index = (2 * row_index_from_top + 1) * row_width + 2 * cell[0] + 1
            result[index] = label

        # Legend at bottom
        for block, label in sorted(labels.items(), key=lambda x: x[1]):
            basis_el, positive = block
            result.append(label)
            result.append(": ")
            if basis_el == (Chord((0,0,1, 1)), Chord((0,1,0,1))) and positive:
                result.append("point")
            else:
                result.append(
                    f"Av{'+' if positive else ''}"
                    f"({', '.join(str(p) for p in basis_el)})"
                )
            result.append("\n")

        if any(not ob.is_single_cell() for ob in self.obstructions):
            result.append("Crossing obstructions:\n")
            for ob in self.obstructions:
                if not ob.is_single_cell():
                    result.append(str(ob))
                    result.append("\n")
        for i, req in enumerate(self.requirements):
            result.append(f"Requirement {i}:\n")
            for r in req:
                result.append(str(r))
                result.append("\n")
        for i, ass in enumerate(self.assumptions):
            result.append(f"Assumption {i}:\n")
            result.append(str(ass))
            result.append("\n")
        if self.assumptions or self.requirements:
            result = result[:-1]

        return "".join(result)'''
        dimensions = "Dimensions: (" + str(self._length) + ", " + str(self._height) + ")\n"

        obs_str = "Obstructions: "
        obs_indent_len = len(obs_str)
        for obstruction in self._obstructions:
            obs_str += str(obstruction)
            obs_str += ",\n" + " " * obs_indent_len
        if len(self._obstructions) == 0:
            obs_str = obs_str[:-1] + "\n"
        else:
            obs_str = obs_str[:-2-obs_indent_len] + "\n"

        reqs_str = "Requirements: "
        reqs_indent_len = len(reqs_str)
        for req_list in self._requirements:
            reqs_str += "{"
            for req in req_list:
                reqs_str += str(req) + "; "
            reqs_str = reqs_str[:-2] + "},\n" + " " * reqs_indent_len
        if len(self._requirements) == 0:
            reqs_str = reqs_str[:-1] + "\n"
        else:
            reqs_str = reqs_str[:-2-reqs_indent_len] + "\n"

        link_str = "Linkages: "
        links_indent_len = len(link_str)
        for linkage in self._linkages:
            link_str += "{"
            for coord in linkage:
                link_str += str(coord) + ", "
            link_str = link_str[:-2] + "},\n" + " " * links_indent_len
        if len(self.linkages) == 0:
            link_str = link_str[:-1]
        else:
            link_str = link_str[:-2-links_indent_len]

        return(dimensions + obs_str + reqs_str + link_str)
    
    @classmethod
    def from_dict(cls, d: dict) -> "Tiling":
        # absolutely no idea if this is right
        """Returns a Tiling object from a dictionary loaded from a JSON
        serialized Tiling object."""
        obstructions = map(GriddedChord.from_dict, d["obstructions"])
        requirements = map(lambda x: map(GriddedChord.from_dict, x), d["requirements"])
        linkages =  map(lambda x: map(tuple, x), d["linkages"])
        assumptions = map(TrackingAssumption.from_dict, d.get("assumptions", []))
        dimensions = d.get("dimensions")
        return cls(
            obstructions=obstructions,
            requirements=requirements,
            linkages=linkages,
            dimensions=dimensions,
            assumptions=assumptions,
        )
