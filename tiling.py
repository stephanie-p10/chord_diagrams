import json
from itertools import chain, filterfalse, product
from typing import (Any, Callable, Dict, FrozenSet, Iterable, Iterator, List, Optional, Set, Tuple)
from collections import Counter, defaultdict

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
)
from .assumptions import TrackingAssumption
from algorithms.map import RowColMap

__all__ = ["Tiling"]

Cell = Tuple[int, int]
ReqList = Tuple[GriddedChord, ...]
CellBasis = Dict[Cell, Tuple[List[Chord], List[Chord]]]
CellFrozenSet = FrozenSet[Cell]
Dimension = Tuple[int, int]
GCTuple = Tuple[GriddedChord, ...]


class Tiling(CombinatorialClass):
    """Tiling class.

    Zero-indexed coordinates/cells from bottom left corner where the (x, y)
    cell is the cell in the x-th column and y-th row.

    Tilings store the obstructions, requirements and linkages but also caches the empty
    cells and the active cells.
    """
    def __init__(self,
        obstructions: Iterable[GriddedChord] = tuple(),
        #point_obstructions: Iterable[Cell] = tuple(),
        requirements: Iterable[Iterable[GriddedChord]] = tuple(), # might need a case for requirement of having one point in a cell...
        #point_requirements: Iterable[Cell] = tuple(),
        linkages: Iterable[Iterable[Cell]] = tuple(),
        assumptions: Iterable[TrackingAssumption] = tuple()):

        super().__init__()
        self._linkages = linkages
        self._obstructions = obstructions
        #self._point_obs = point_obstructions
        self._requirements = requirements
        #self._point_reqs = point_requirements
        self._assumptions = assumptions
        #self.cells = []

        self._cached_properties = {}
        self._cached_properties["forward_map"] = RowColMap.identity((0, 0))
        self._cached_properties["backward_map"] = RowColMap.identity((0, 0))

        # currently defaults cells with no requirements or obsturctions to be empty
        # note: should this be defaulted to allowing other chords?
        self._prepare_properties()

    # also changed how acitve_cells and empty_cells are computed, no longer add point obs based on empty cells
    def _prepare_properties(self) -> None:
        """
        Compute _active_cells, _empty_cells, _dimensions, and store them
        """
        active_cells = union_reduce(
            set(ob.pos) for ob in self._obstructions if not ob.is_point()
        )
        active_cells.update(
            *(union_reduce(set(comp.pos) for comp in req) for req in self._requirements)
        )

        max_row = 0
        max_col = 0
        for cell in active_cells:
            max_col = max(max_col, cell[0])
            max_row = max(max_row, cell[1])
        dimensions = (max_col + 1, max_row + 1)

        empty_cells = tuple(
            cell
            for cell in product(range(dimensions[0]), range(dimensions[1]))
            if cell not in active_cells
        )

        # If the first obstruction is the empty perm, we shouldn't do any of this.
        if len(self._obstructions) > 0 and len(self._obstructions[0]) > 0:
            # We can assume that self._obstructions is sorted at this point, so to
            #   extract the point obstructions, we just pass though them until we've
            #   found the last one, then we slice the list there.
            self._obstructions = sorted(self._obstructions) # not assuming, fix later sToDo
            index = 0
            for ob in self._obstructions:
                if len(ob) > 1:
                    break
                index += 1  # Now the last point obstruction is at index [index-1]
            non_point_obstructions = self._obstructions[index:]

            new_point_obstructions = tuple(
                GriddedChord(Chord(0,), (cell,)) for cell in empty_cells
            )
            self._obstructions = new_point_obstructions + non_point_obstructions

        self._cached_properties["active_cells"] = frozenset(active_cells)
        self._cached_properties["empty_cells"] = frozenset(empty_cells)
        self._cached_properties["dimensions"] = dimensions

    def _minimize_mapping(self) -> RowColMap:
        """
        Returns a pair of dictionaries, that map rows/columns to an
        equivalent set of rows/columns where empty ones have been removed.
        Also returns a boolean describing whether this mapping is the identity
        mapping which saves some later computation.
        """
        active_cells = self.active_cells
        assert active_cells
        col_set = set(c[0] for c in active_cells)
        row_set = set(c[1] for c in active_cells)
        col_list, row_list = sorted(col_set), sorted(row_set)
        identity = (self.dimensions[0] == len(col_list)) and (
            self.dimensions[1] == len(row_list)
        )
        col_mapping = {x: actual for actual, x in enumerate(col_list)}
        row_mapping = {y: actual for actual, y in enumerate(row_list)}
        return RowColMap(row_map=row_mapping, col_map=col_mapping, is_identity=identity)

    def _remove_empty_rows_and_cols(self) -> None:
        """Remove empty rows and columns."""
        # Produce the mapping between the two tilings
        if not self.active_cells:
            assert GriddedChord.empty_chord() not in self.obstructions
            self._cached_properties["forward_map"] = RowColMap.identity((0, 0))
            self._obstructions = (GriddedChord.single_cell((0,), (0, 0)),)
            self._requirements = tuple()
            self._assumptions = tuple()
            self._cached_properties["dimensions"] = (1, 1)
            return
        forward_map = self._minimize_mapping()
        self._cached_properties["forward_map"] = forward_map
        # We still may need to remove point obstructions if the empty row or col
        # was on the end so we do it outside the next if statement.
        self._obstructions = tuple(
            forward_map.map_gc(ob)
            for ob in self.obstructions
            if not ob.is_point() or forward_map.is_mappable_gc(ob)
        )

        if not forward_map.is_identity():
            self._requirements = tuple(
                tuple(forward_map.map_gc(req) for req in reqlist)
                for reqlist in self._requirements
            )
            self._assumptions = tuple(
                sorted(
                    forward_map.map_assumption(assumption)
                    for assumption in self._assumptions
                )
            )
            self._linkages = tuple(
                tuple(forward_map.map_cell(cell) for cell in linkage if forward_map.is_mappable_cell(cell))
                for linkage in self._linkages
            )
            self._cached_properties["active_cells"] = frozenset(
                forward_map.map_cell(cell)
                for cell in self._cached_properties["active_cells"]
                if forward_map.is_mappable_cell(cell)
            )
            self._cached_properties["empty_cells"] = frozenset(
                forward_map.map_cell(cell)
                for cell in self._cached_properties["empty_cells"]
                if forward_map.is_mappable_cell(cell)
            )
            self._cached_properties["dimensions"] = (
                forward_map.max_col() + 1,
                forward_map.max_row() + 1,
            )

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
    
    @property
    def dimensions(self):
        try:
            dims = (self._cached_properties["dimensions"])
            return dims
        except KeyError:
            self._prepare_properties()
            return self._cached_properties["dimensions"]

    @property
    def forward_map(self) -> RowColMap:
        try:
            return self._cached_properties["forward_map"]
        except KeyError:
            self._remove_empty_rows_and_cols()
            return self._cached_properties["forward_map"]

    @property
    def backward_map(self) -> RowColMap:
        try:
            return self._cached_properties["backward_map"]
        except KeyError:
            backward_map = self.forward_map.reverse()
            self._cached_properties["backward_map"] = backward_map
            return backward_map
    
    @property
    def active_cells(self) -> CellFrozenSet:
        """
        Returns a set of all cells that do not contain a point obstruction,
        i.e., not empty.
        """
        try:
            return self._cached_properties["active_cells"]
        except KeyError:
            self._prepare_properties()
            return self._cached_properties["active_cells"]
        
    # I don't think we care about linkages; if a linkage has a row of empty cells, then we can still reduce.
    def _remove_empty_rows_and_cols(self) -> None:
        """Remove empty rows and columns."""
        # Produce the mapping between the two tilings
        if not self.active_cells:
            assert GriddedChord.empty_chord() not in self.obstructions
            self._cached_properties["forward_map"] = RowColMap.identity((0, 0))
            self._obstructions = (GriddedChord.single_cell(Chord(0,), (0, 0)),)
            self._requirements = tuple()
            self._assumptions = tuple()
            self._cached_properties["dimensions"] = (1, 1)
            return
        forward_map = self._minimize_mapping()
        self._cached_properties["forward_map"] = forward_map
        # We still may need to remove point obstructions if the empty row or col
        # was on the end so we do it outside the next if statement.
        self._obstructions = tuple(
            forward_map.map_gp(ob)
            for ob in self.obstructions
            if not ob.is_point_perm() or forward_map.is_mappable_gp(ob)
        )

        if not forward_map.is_identity():
            self._requirements = tuple(
                tuple(forward_map.map_gp(req) for req in reqlist)
                for reqlist in self._requirements
            )
            self._assumptions = tuple(
                sorted(
                    forward_map.map_assumption(assumption)
                    for assumption in self._assumptions
                )
            )
            self._cached_properties["active_cells"] = frozenset(
                forward_map.map_cell(cell)
                for cell in self._cached_properties["active_cells"]
                if forward_map.is_mappable_cell(cell)
            )
            self._cached_properties["empty_cells"] = frozenset(
                forward_map.map_cell(cell)
                for cell in self._cached_properties["empty_cells"]
                if forward_map.is_mappable_cell(cell)
            )
            self._cached_properties["dimensions"] = (
                forward_map.max_col() + 1,
                forward_map.max_row() + 1,
            )

    def cells_in_row(self, row: int) -> CellFrozenSet:
        """Return all active cells in row."""
        return frozenset((x, y) for (x, y) in self.active_cells if y == row)

    def cells_in_col(self, col: int) -> CellFrozenSet:
        """Return all active cells in column."""
        return frozenset((x, y) for (x, y) in self.active_cells if x == col)
    
    def get_assumption_parameter(self, assumption: TrackingAssumption) -> str:
        """
        Return the variable associated with the given assumption.

        Raise ValueError if the assumptions is not on the tiling.
        """
        try:
            idx = tuple(self._assumptions).index(assumption)
        except ValueError as e:
            raise ValueError(
                f"following assumption not on tiling: '{assumption}'"
            ) from e
        return f"k_{idx}"

    # sTODO this is currently incorrect, but ok for non crossing. the atom should be a set containing the smallest thing that can be made.
    def is_atom(self):
        """Return True if the Tiling is a single gridded chord."""
        single_size_1 = False
        size_1_in_self = []
        cells = self._cached_properties["active cells"]
        for gc in GriddedChord.all_grids(Chord((1,1)), cells):
            if self.contains(gc):
                size_1_in_self.append(gc)

        single_size_1 = (len(size_1_in_self) == 1)
        
        size_2_chords = [Chord((0,0,1,1)), Chord((0,1,0,1)), Chord((0,1,1,0))]
        avoids_all_2s = True
        for chord in size_2_chords:
            for gc in GriddedChord.all_grids(chord, cells):
                if self.contains(gc):
                    avoids_all_2s = False
                    break
        return single_size_1 and avoids_all_2s

    # currently this is almost the same as is_empty, should probably also be optimized.
    def minimum_size_of_object(self):
        """Return the size of the smallest object in the class."""
        # finds sum of the maximum length requirements in each requirements list
        sum_max_reqs = 0
        for req_list in self.requirements:
            if len(req_list) == 0:
                return 0
            req_lengths = [len(req) for req in req_list]
            sum_max_reqs += max(req_lengths)

        if sum_max_reqs == 0: # then there are no requirements, but we still need a chord of size one
            sum_max_reqs += 1
        
        # proved maximum size of smallest chord that can be gridded:
        max_len = sum_max_reqs * 2 - 1

        all_chords = []
        for num_chords in range(1, max_len + 1):
            all_chords += list(Chord.of_length(num_chords))

        # product of length and width to get valid cells can probaby be much improved.
        cells = self._cached_properties["active cells"]
        for chord in all_chords:
            for gc in GriddedChord.all_grids(chord, cells):
                if self.contains(gc):
                    return len(gc)
        
        return 0
    
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
       
        all_chords = []
        for num_chords in range(1, max_len + 1):
            all_chords += list(Chord.of_length(num_chords))

        # product of length and width to get valid cells can probaby be much improved.
        cells = self._cached_properties["active cells"]
        for chord in all_chords:
            for gc in GriddedChord.all_grids(chord, cells):
                if self.contains(gc):
                    return False

        return True

    def contains(self, gc: GriddedChord) -> bool:
        has_reqs = all(gc.contains(*req) for req in self._requirements)
        avoids_ob = not any(gc.contains(ob) for ob in self._obstructions)
        links_connected = all(gc.is_connected(cells) for cells in self._linkages)
        return has_reqs and avoids_ob and links_connected
    
    def add_list_requirement(self, req_list: Iterable[GriddedChord]) -> "Tiling":
        """
        Return a new tiling with the requirement list added.
        """
        new_req = tuple(sorted(req_list))
        return Tiling(
            self._obstructions,
            sorted(self._requirements + (new_req,)),
            self._linkages,
            self._assumptions,
        )

    def add_obstructions(self, gps: Iterable[GriddedChord]) -> "Tiling":
        """Returns a new tiling with the obstructions added."""
        new_obs = tuple(gps)
        return Tiling(
            tuple(self._cached_properties["dimensions"][0], self._cached_properties["dimensions"][1]),
            sorted(self._obstructions + new_obs),
            self._requirements,
            self._linkages,
            self._assumptions
        )
    
    def only_cell_in_col(self, cell: Cell) -> bool:
        """Checks if the cell is the only active cell in the column."""
        return sum(1 for (x, y) in self.active_cells if x == cell[0]) == 1

    def only_cell_in_row(self, cell: Cell) -> bool:
        """Checks if the cell is the only active cell in the row (besides other end of chord)."""
        return sum(1 for (x, y) in self.active_cells if y == cell[1]) == 2
    
    #sTODO this is wrong! currently this works for when obstructions are simplified automatically.
    @property
    def point_cells(self) -> CellFrozenSet:
        try:
            return self._cached_properties["point_cells"]
        except KeyError:
            # finds all cells obstructing length one and two chords fully conatined within them
            local_length_lt2_obcells = Counter(
                ob.pos[0]
                for ob in self._obstructions
                if ob.is_localized() and (ob._patt == (0, 0) or ob._patt == (0, 1))
            )
            # finds cells that must only have a single point
            point_cells = frozenset(
                cell for cell in self.positive_cells if local_length_lt2_obcells[cell] == 2
            )
            self._cached_properties["point_cells"] = point_cells
            return point_cells
        
    @property
    def positive_cells(self) -> CellFrozenSet:
        """Cells that must have something in them"""
        try:
            return self._cached_properties["positive_cells"]
        except KeyError:
            positive_cells = frozenset(
                union_reduce(
                    intersection_reduce(req.pos for req in reqs)
                    for reqs in self._requirements
                )
            )
            self._cached_properties["positive_cells"] = positive_cells
            return positive_cells
    
    @property
    def active_cells(self) -> CellFrozenSet:
        """
        Returns a set of all cells that do not contain a point obstruction,
        i.e., not empty.
        """
        try:
            return self._cached_properties["active_cells"]
        except KeyError:
            self._prepare_properties()
            return self._cached_properties["active_cells"]

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
        dimensions = "Dimensions: (" + str(self._cached_properties["dimensions"][0]) + ", " + str(self._cached_properties["dimensions"][1]) + ")\n"

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
    
    def to_jsonable(self):
        output: dict = super().to_jsonable()
        output["obstructions"] = [gc.to_jsonable() for gc in self.obstructions]
        output["requirements"] = [[gc.to_jsonable() for gc in req] for req in self.requirements]
        output["linkages"] = self._linkages
        output["dimensions"] = tuple([self._cached_properties["dimensions"][0], self._cached_properties["dimensions"][1]])
        output["assumptions"] = [assump.to_jsonable() for assump in self.assumptions]
        return output

    @classmethod
    def from_dict(cls, d: dict) -> "Tiling":
        # reasonably sure this is correct, not sure about formatting
        """Returns a Tiling object from a dictionary loaded from a JSON
        serialized Tiling object."""
        obstructions = tuple(map(GriddedChord.from_dict, d["obstructions"]))
        requirements = tuple(map(lambda x: tuple(map(GriddedChord.from_dict, x)), d["requirements"]))
        linkages =  tuple(map(lambda x: tuple(map(tuple, x)), d["linkages"]))
        assumptions = tuple(map(TrackingAssumption.from_dict, d.get("assumptions", [])))
        dimensions = d.get("dimensions")
        return cls(
            obstructions=tuple(obstructions),
            requirements=tuple(requirements),
            linkages=tuple(linkages),
            dimensions=dimensions,
            assumptions=tuple(assumptions),
        )
