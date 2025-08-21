import json
from itertools import chain, filterfalse, product, combinations
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
from assumptions import TrackingAssumption
from algorithms.map import RowColMap
from algorithms.simplify import SimplifyObstructionsAndRequirements
from algorithms.expansion import Expansion

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
        requirements: Iterable[Iterable[GriddedChord]] = tuple(), 
        linkages: Iterable[Iterable[Cell]] = tuple(),
        assumptions: Iterable[TrackingAssumption] = tuple(),
        simplify: bool = True,
        remove_empty_rows_and_cols: bool = True,
        derive_empty: bool = True,
        expand: bool = True):

        super().__init__()
        self._linkages = tuple(linkages)
        self._obstructions = tuple(obstructions)
        self._requirements = tuple(requirements)
        self._assumptions = tuple(assumptions)
        self._cached_properties = {}

        if "dimensions" not in self._cached_properties:
            self._compute_dimensions()

        #print("computed dimensions")

        if simplify:
            self._simplify()

        #print("computed simplify 1")
        
        # currently defaults cells with no requirements or obsturctions to be empty
        # note: should this be defaulted to allowing other chords?
        if "empty_cells" not in self._cached_properties:
            self._prepare_properties(derive_empty)

        #print("computed prepare properties")

        if expand:
            self._expand()

        #print("computed expand")

        # for simplification of computation, point obstructions to mark empty cells are added after expand.
        # this should not affect the switch to only describe a tiling with chord diagrams, since the only time
        #   the new point obstructions are used is when dealing with parts of the tiling that are empty. Most 
        #   computations are done only on the active cells.
        #print("empty cells", self._cached_properties["empty_cells"])
        #self._add_point_obs(self._cached_properties["empty_cells"])

        #print("computed adding point obs")
            
        if simplify:
            self._simplify()

        #print("computed simplify 2")

        if remove_empty_rows_and_cols:
            self._remove_empty_rows_and_cols()

        #print("computed remove empty rows and cols")
        #print()

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
    def dimensions(self) -> Tuple[int, int]:
        try:
            dims = (self._cached_properties["dimensions"])
            return dims
        except KeyError:
            self._compute_dimensions()
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
        
    #s TODO this is wrong! currently this works for allowing non chord patterns
    @property
    def point_cells(self) -> CellFrozenSet:
        try:
            return self._cached_properties["point_cells"]
        except KeyError:
            # finds all cells obstructing length one and two chords fully contained within them
            local_length_lt2_obcells = Counter(
                ob.pos[0]
                for ob in self._obstructions
                if ob.is_localized() and (ob._patt == (0, 0) or ob._patt == (0, 1) or ob._patt == (1,0))
            )
            # finds cells that must only have a single point
            point_cells = frozenset(
                cell for cell in self.positive_cells if local_length_lt2_obcells[cell] == 3
            )
            self._cached_properties["point_cells"] = point_cells
            return point_cells
        
    @property
    def chord_cells(self) -> CellFrozenSet:
        try:
            return self._cached_properties["chord_cells"]
        except KeyError:
            local_len_2_obcells = Counter(
                ob.pos[0]
                for ob in self._obstructions
                if ob.is_localized() and (ob._patt == (0,1) or ob._patt == (1,0))
            )
            chord_cells = frozenset(cell for cell in self.positive_cells if local_len_2_obcells[cell] == 2)
            self._cached_properties["chord_cells"] = chord_cells
            return chord_cells
        
    def _compute_dimensions(self) -> None:
        max_x = 0
        max_y = 0

        for ob in self._obstructions:
            for x, y in ob.pos:
                if x > max_x:
                    max_x = x
                if y > max_y:
                    max_y = y

        for reqlist in self._requirements:
            for req in reqlist:
                for x, y in req.pos:
                    if x > max_x:
                        max_x = x
                    if y > max_y:
                        max_y = y

        for link in self._linkages:
            for (x, y) in link:
                if x > max_x:
                    max_x = x
                if y > max_y:
                    max_y = y

        self._cached_properties["dimensions"] = (max_x + 1, max_y + 1)

        return 

    
    # also changed how acitve_cells and empty_cells are computed, no longer add point obs based on empty cells
    def _prepare_properties(self, derive_empty: bool = True) -> None:
        """
        Compute active_cells and empty_cells, and store them in cached_properties
        """
        potential_active_cells = []
        max_x = self.dimensions[0]
        max_y = self.dimensions[1]
        all_cells = list(product(range(max_x), range(max_y)))

        # If we are derviving the tiling empty, a cell can only be active if it is found in an ob or req
        if derive_empty:
            cell_set = set()
            cell_set = union_reduce(
                set(ob.pos) for ob in self._obstructions if not ob.is_point()
            )
            cell_set.update(
                *(union_reduce(set(comp.pos) for comp in req) for req in self._requirements)
            )
            potential_active_cells = list(cell_set)
        # If we were not assuming everything non used cell is empty, any cell could be an active cell
        else:
            potential_active_cells = all_cells

        # A cell is empty if there are no chord diagrams that can be gridded on the tiling with a point in the cell
        # If there is a chord diagram with a point in the cell, there will be one with at most size max_size_to_check
        max_size_to_check = self.maximum_length_of_minimum_gridded_chord() + 1
        #print(max_size_to_check)
        all_chords = []
        for num_chords in range(1, max_size_to_check + 1):
            all_chords += list(Chord.of_length(num_chords))

        # Checks what cells are used in any gridded chord diagram up to max_size_to_check
        cells_used = set()
        for chord in all_chords:
            for gc in GriddedChord.all_grids(chord, potential_active_cells):
                if self.contains(gc):
                    cells_used.update(gc.pos)

        # calculates empty cells as complement of active cells
        empty_cells = tuple(
            cell
            for cell in all_cells
            if cell not in cells_used
        )

        self._cached_properties["active_cells"] = frozenset(cells_used)
        self._cached_properties["empty_cells"] = frozenset(empty_cells)

    def _add_point_obs(self, cells: Tuple[Cell]):
        #print("add_point_obs")
        new_obs = list(self._obstructions)
        #print("new_obs", new_obs)
        #print("cells", cells)
        for cell in cells:
            new_obs.append(GriddedChord.single_cell(Chord((0,)), cell))
            #print(cell)

        self._obstructions = tuple(new_obs)
  
    def _expand(self) -> None:
        expansion_class = Expansion(self._obstructions, self._requirements, self.dimensions, self.active_cells)
        #print("printing obs in expand", self._obstructions)
        #print("printing reqs in expand", self._requirements)
        expansion_class.expand_obstructions()
        expansion_class.expand_requirements()
        self._obstructions = tuple(expansion_class._obstructions)
        self._requirements = tuple(expansion_class._requirements)
        #print("printing obs in expand after algo", self._obstructions)
        #print("printing reqs in expand after algo", self._requirements)

    def _simplify(self) -> None:
        simplify_algo = SimplifyObstructionsAndRequirements(self.obstructions, self.requirements, self._cached_properties["dimensions"])
        simplify_algo.simplify()
        self._obstructions = simplify_algo.obstructions
        self._requirements = simplify_algo.requirements

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

    # I don't think we care about linkages; if a linkage has a row of empty cells, then we can still reduce.
    def _remove_empty_rows_and_cols(self) -> None:
        """Remove empty rows and columns."""
        # Produce the mapping between the two tilings
        if not self.active_cells:
            assert GriddedChord.empty_chord() not in self.obstructions
            self._cached_properties["forward_map"] = RowColMap.identity((0, 0))
            self._obstructions = (GriddedChord.single_cell(Chord((0,0)), (0, 0)),)
            self._requirements = tuple()
            self._assumptions = tuple()
            self._cached_properties["dimensions"] = (1, 1)
            return
        forward_map = self._minimize_mapping()
        self._cached_properties["forward_map"] = forward_map
        # We still may need to remove point obstructions if the empty row or col
        # was on the end so we do it outside the next if statement.
        #print(forward_map)
        self._obstructions = tuple(
            forward_map.map_gc(ob)
            for ob in self.obstructions
            if not ob.is_point() and forward_map.is_mappable_gc(ob)
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

    # proved bound has -1 at the end, but I'm not sure how to account for reqs that are not valid chord diagrams
    def maximum_length_of_minimum_gridded_chord(self) -> int:
        """Returns the maximum length of the minimum gridded chord diagram that
        can be gridded on the tiling.
        """
        if not self.requirements:
            #print(self.requirements)
            return 1
        
        sum_largest_req_lengths = 0
        for reqs in self.requirements:
            if len(reqs) != 0:
                req_lengths = [max(req.patt) + 1 for req in reqs if len(req.patt) > 0]
                sum_largest_req_lengths += max(req_lengths)
        
        return 2 * sum_largest_req_lengths - 1

    def sub_tiling(
        self,
        cells: Iterable[Cell],
        factors: bool = False,
        add_assumptions: Iterable[TrackingAssumption] = tuple(),
    ) -> "Tiling":
        """Return the tiling using only the obstructions and requirements
        completely contained in the given cells. If factors is set to True,
        then it assumes that the first cells confirms if a gridded perm uses only
        the cells."""
        obstructions = tuple(
            ob
            for ob in self.obstructions
            if (factors and ob.pos[0] in cells) or all(c in cells for c in ob.pos)
        )
        requirements = tuple( # tuple -> Tiling.sort_requirements when implemented sToDo
            req
            for req in self.requirements
            if (factors and req[0].pos[0] in cells)
            or all(c in cells for c in chain.from_iterable(r.pos for r in req))
        )
        assumptions = tuple(
            assump.__class__(
                gc
                for gc in assump.gcs
                if (factors and gc.pos[0] in cells) or all(c in cells for c in gc.pos)
            )
            for assump in self.assumptions
        ) + tuple(add_assumptions)
        linkages = tuple(
            link
            for link in self.linkages
            if (factors and link[0] in cells)
            or all(cell in cells for cell in link)
        )
        
        # TODO: check sum/skew assumptions (this comment was inherited)
        return self.__class__(
            obstructions,
            requirements,
            linkages,
            tuple(sorted(set(assump for assump in assumptions if assump.gcs))),
            simplify=False,
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

    # sToDo this is currently incorrect, but ok for non crossing. the atom should be a set containing the smallest thing that can be made.
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
    
    # sToDo this is inefficient and in permutations there is a class that does this efficiently
    # gridded_perms_of_length -> all_chords_on_tiling
    def all_chords_on_tiling(self, size: int = 0, use_non_chord_patts: bool = False) -> list[GriddedChord]:
        """Returns all patterns with up to size points that can be gridded on the tiling"""
        all_chords = []
        for num_chords in range(1, size//2 + 1):
            all_chords += list(Chord.of_length(num_chords))
        
        chords_on_tiling = set() # this is a set so the non chord diagram pattern code works nicely
        cells = self.active_cells
        #print(all_chords)
            
        for chord in all_chords:
            #print(list(GriddedChord.all_grids(chord, cells)))
            #print(list(filter(self.contains, list(GriddedChord.all_grids(chord, cells)))))
            #print([self.contains(gc) for gc in list(GriddedChord.all_grids(chord, cells))])
            chords_on_tiling.update(filter(self.contains, GriddedChord.all_grids(chord, cells)))
        
        #print(chords_on_tiling)

        if use_non_chord_patts:
            bigger_chords = []
            for num_chords in range(size//2 + 1, size + 1):
                bigger_chords += list(map(lambda chord: chord.get_pattern(), Chord.of_length(num_chords)))
                
            for chord_size in range(1, size + 1):
                for chord in bigger_chords:
                    #print(list(combinations(chord, chord_size)))
                    for subchord in list(combinations(chord, chord_size)):
                        subchord = Chord.to_standard(subchord, True)
                        chords_on_tiling.update(filter(self.contains, GriddedChord.all_grids(subchord, cells)))
                        #if chord_size == 2:
                            #print(list(GriddedChord.all_grids(subchord, cells)))
                            #print(list(filter(self.contains, GriddedChord.all_grids(subchord, cells))))
                    #if chord_size == 2:
                        #print()
                #print(chord_size)
                    
        return list(chords_on_tiling)

    def is_empty(self) -> bool:
        """Checks if the tiling is empty.

        Tiling is empty if it has been inferred to be contradictory due to
        contradicting requirements and obstructions or no gridded chord
        can be gridded on the tiling.
        """
        # if any obstruction is empty, no chords can be gridded since all grids contain the empty grid.
        if any(ob.is_empty() for ob in self.obstructions):
            #print("contains empty obstruction")
            return True
        
        # if any req list is empty or has all empty grids, only the empty grid can be gridded. 
        if any(all(req.is_empty() for req in req_list) for req_list in self.requirements) and len(self.requirements) >= 1:
            #print("contains empty requirements list")
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
        cells = self._cached_properties["active_cells"]
        for chord in all_chords:
            for gc in GriddedChord.all_grids(chord, cells):
                #print(gc)
                if self.contains(gc):
                    #print("found contains", gc)
                    return False

        return True

    def contains(self, gc: GriddedChord) -> bool:
        has_reqs = all(gc.contains(*req) for req in self._requirements)
        avoids_ob = not any(gc.contains(ob) for ob in self._obstructions)
        links_connected = all(gc.is_connected(cells) for cells in self._linkages)
        #print(has_reqs, avoids_ob, links_connected)
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

    def add_requirement(self, patt: Chord, pos: Iterable[Cell]) -> "Tiling":
        """Returns a new tiling with the requirement of the pattern
        patt with position pos."""
        new_req_list = (GriddedChord(patt, pos),)
        return self.add_list_requirement(new_req_list)

    def add_obstructions(self, gcs: Iterable[GriddedChord], simplify: bool = False, expand: bool = False) -> "Tiling":
        """Returns a new tiling with the obstructions added."""
        new_obs = tuple(gcs)
        #print(sorted(self._obstructions + new_obs))
        all_obs = sorted(self._obstructions + new_obs)
        return Tiling(
            all_obs,
            self._requirements,
            self._linkages,
            self._assumptions,
            simplify= simplify,
            expand= expand
        )
    
    def add_assumptions(
        self, assumptions: Iterable[TrackingAssumption], clean: bool = True
    ) -> "Tiling":
        """Returns a new tiling with the added assumptions."""
        tiling = Tiling(
            self._obstructions,
            self._requirements,
            self._linkages,
            self._assumptions + tuple(assumptions),
            remove_empty_rows_and_cols=False,
            derive_empty=False,
            simplify=False,
            sorted_input=True,
        )
        #if clean: ## fix when implemented
        #    tiling.clean_assumptions()
        return tiling
    
    def only_cell_in_col(self, cell: Cell) -> bool:
        """Checks if the cell is the only active cell in the column."""
        return sum(1 for (x, y) in self.active_cells if x == cell[0]) == 1

    def only_cell_in_row(self, cell: Cell) -> bool:
        """Checks if the cell is the only active cell in the row"""
        return sum(1 for (x, y) in self.active_cells if y == cell[1]) == 1
    
    def obs_in_cell(self, cell: Cell) -> tuple[int]:
        """Returns the indices obstructions in cell"""

    def reqs_in_cell(self, cell: Cell) -> tuple[int]:
        """Returns the indices of requirments lists where at least one requirment in the list has a point in cell"""

    def __hash__(self) -> int:
        return (
            hash(self._requirements)
            ^ hash(self._obstructions)
            ^ hash(self._assumptions)
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Tiling):
            return False
        #print(self.obstructions == other.obstructions, 
        #      self.requirements == other.requirements, 
        #      self.linkages == other.linkages, 
        #      self.assumptions == other.assumptions)
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
        return cls(
            obstructions=tuple(obstructions),
            requirements=tuple(requirements),
            linkages=tuple(linkages),
            assumptions=tuple(assumptions),
        )
