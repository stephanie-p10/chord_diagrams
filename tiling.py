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
        #fixed_rows: Iterable[int] = tuple(),
        #fixed_cols: Iterable[int] = tuple(),
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
        self._derive_empty = derive_empty
        #self._fixed_rows = fixed_rows
        #self._fixed_cols = fixed_cols

        if "dimensions" not in self._cached_properties:
            self._compute_dimensions()

        if simplify:
            self._simplify()

        # currently defaults cells with no requirements or obsturctions to be empty
        # note: should this be defaulted to allowing other chords?
        if "empty_cells" not in self._cached_properties and derive_empty:
            #print("used empty cells")
            self._prepare_properties(derive_empty)

        if expand:
            self._expand()

        #print("computed expand")

        # for simplification of computation, point obstructions to mark empty cells are added after expand.
        # this should not affect the switch to only describe a tiling with chord diagrams, since the only time
        #   the new point obstructions are used is when dealing with parts of the tiling that are empty. Most 
        #   computations are done only on the active cells.
        if derive_empty:
            self._add_point_obs(self._cached_properties["empty_cells"])

        #print("computed adding point obs")
            
        # this will take out all obs that have ends in empty cells
        if simplify:
            self._simplify()

        if remove_empty_rows_and_cols:
            self._remove_empty_rows_and_cols()


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
        
    # sTODO tests. I changed this so it works only for chords that are described in full chord diagrams
    @property
    def point_cells(self) -> CellFrozenSet:
        """ Finds cells that cannot have more than one point in them
        """
        try:
            return self._cached_properties["point_cells"]
        except KeyError:
            potential_cells = list(self.active_cells)

            # generates all chords on the tiling of length two
            size_two_chords = self.chords_of_length(2)

            # if a cell has more than two points in it, it must contain the pattern 00, 01, or 10
            # this means that the cells that can have one point are the cells that do not contain 00, 01 or 10
            # in other words, they do not contain 00, or a size two chord pattern with at least two points in that cell
            bad_cells = []

            # finding which cells have 00
            for cell in potential_cells:
                if self.contains(GriddedChord(Chord((0, 0)), (cell, cell))):
                    bad_cells.append(cell)

            # finding the cells that contain two points of a size two chord
            for gc in size_two_chords:
                pos_list_counter = Counter(gc.pos)
                for cell in pos_list_counter.keys():
                    if pos_list_counter[cell] > 1:
                        bad_cells.append(cell)

            required_cells = []
            for reqlist in self.requirements:
                cells_in_all_reqs = intersection_reduce(req.pos for req in reqlist)
                required_cells += cells_in_all_reqs

            required_cells_set = set(required_cells)

            # filtering out all the cells that contain more than one point
            point_cells = []
            for cell in potential_cells:
                if cell not in bad_cells and cell in required_cells_set:
                    point_cells.append(cell)

            self._cached_properties["point_cells"] = point_cells
            return point_cells
        
    @property
    # sTODO TESTS!!
    def chord_cells(self) -> CellFrozenSet:
        """Not tested, but should retrun cells that have exactly one chord in them"""
        try:
            return self._cached_properties["chord_cells"]
        except KeyError:
            # chord cells are cells that contain exactly one chord (both ends)
            # this means they contain 00, and avoid 01, 10
            possible_00_chords = [GriddedChord(Chord((0, 0)), (cell, cell)) for cell in self.active_cells]
            cells_with_00 = []
            
            # check which cells must contain a 00 chord diagram pattern
            for reqlist in self.requirements:
                if len(reqlist) > 1:
                    for single_chord in possible_00_chords:
                        if all(req.contains(single_chord) for req in reqlist):
                            cell = single_chord.pos[0]
                            cells_with_00.append(cell)

            exp_class = Expansion(self.obstructions, self.requirements, self._cached_properties["dimensions"])
            good_cells = cells_with_00
            # of the cells that contain 00, checks which ones do not contain a 01 and 10 pattern
            for cell in cells_with_00:
                # this is a terrible variable name
                _01_10_chords_in_cell = exp_class.expand_gridded_chords([GriddedChord(Chord((0, 1)), (cell, cell)),
                                                                         GriddedChord(Chord((1, 0)), (cell, cell))],
                                                                         self.obstructions)
                if len(_01_10_chords_in_cell) != 0:
                    good_cells.remove(cell)
            return good_cells
        
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

    def _prepare_properties(self, derive_empty = True) -> None:
        """
        Compute active_cells and empty_cells, and store them in cached_properties
        """
        ''' OLD VERY SLOW CODE (checks all edge cases)
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
        )'''

        # Fast method of calculating active cells, assumes the user did not do anything "silly"
        # adds all cells that an obstuction larger than a single point uses
        active_cells = union_reduce(
            set(ob.pos) for ob in self.obstructions if len(ob.patt) > 1
        )
        # adds all cells that are used in a requirement to the active cells.
        active_cells.update(
            *(union_reduce(set(comp.pos) for comp in req) for req in self.requirements)
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

        if self._derive_empty:
            self._cached_properties["active_cells"] = frozenset(active_cells)
            self._cached_properties["empty_cells"] = frozenset(empty_cells)
            self._cached_properties["dimensions"] = dimensions
        else:
            self._cached_properties["active_cells"] = list(product(range(max_col + 1), range(max_row + 1)))
            self._cached_properties["empty_cells"] = frozenset()
            self._cached_properties["dimensions"] = dimensions

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
        expansion_class.expand_obstructions()
        expansion_class.expand_requirements()
        self._obstructions = tuple(expansion_class._obstructions)
        self._requirements = tuple(expansion_class._requirements)

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

    def objects_of_size(self, size):
        print("called objects of size from us")
        chords = Chord.of_length(size)

        grids = []
        for chord in chords:
            for grid in GriddedChord.all_grids(chord):
                if self.contains(grid):
                    grids.append(grid)

        return grids

    # sToDo this is currently incorrect, but ok for non crossing. the atom should be a set containing the smallest thing that can be made.
    def is_atom(self):
        """Return True if the Tiling is a single gridded chord."""
        
        # if any req list is empty or has all empty grids, only the empty grid can be gridded. 
        if any(all(req.is_empty() for req in req_list) and 
               len(req_list) != 0 for req_list in self.requirements) and len(self.requirements) >= 1:
            #print("contains empty requirements list")
            return True

        single_size_1 = False
        size_1_in_self = []
        try:
            cells = self._cached_properties["active_cells"]
        except(KeyError):
            self._prepare_properties()
            cells = self._cached_properties["active_cells"]

        for gc in GriddedChord.all_grids(Chord((1,1)), cells):
            if self.contains(gc):
                size_1_in_self.append(gc)

        if len(size_1_in_self) == 0 and self.contains(GriddedChord(Chord(()), ())):
            return True

        single_size_1 = (len(size_1_in_self) == 1)
        
        size_2_chords = [Chord((0,0,1,1)), Chord((0,1,0,1)), Chord((0,1,1,0))]
        avoids_all_2s = True
        for chord in size_2_chords:
            for gc in GriddedChord.all_grids(chord, cells):
                if self.contains(gc):
                    avoids_all_2s = False
                    break
        return single_size_1 and avoids_all_2s and not self.contains(GriddedChord(Chord(()), ()))

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

        # proved maximum size of smallest chord that can be gridded:
        max_len = sum_max_reqs * 2 - 1

        if sum_max_reqs == -1: # then there are no requirements, but we still need a chord of size one
            sum_max_reqs = 0

        all_chords = []
        for num_chords in range(0, max_len + 1):
            all_chords += list(Chord.of_length(num_chords))

        # product of length and width to get valid cells can probaby be much improved.
        try:
            cells = self._cached_properties["active_cells"]
        except KeyError:
            self._prepare_properties()
            cells = self._cached_properties["active_cells"]
        
        for chord in all_chords:
            for gc in GriddedChord.all_grids(chord, cells):
                if self.contains(gc):
                    return len(gc)
        
        return 0 # this maybe should be -1?
    
    # s TODO fix all_chords_on_tiling to use this
    def chords_of_length(self, length: int, use_non_chord_patts: bool = False) -> list[GriddedChord]:
        """Returns all chord diagrams on the tiling of length length"""
        all_chords = list(Chord.of_length(length))
        cells = self.active_cells
        chords_on_tiling = set()
        for chord in all_chords:
            chords_on_tiling.update(filter(self.contains, GriddedChord.all_grids(chord, cells)))

        return list(chords_on_tiling)

    
    # sToDo this is inefficient and in permutations there is a class that does this efficiently
    # gridded_perms_of_length -> all_chords_on_tiling
    # s TODO check where this is used, see if the empty chord should be included?
    def all_chords_on_tiling(self, size: int = 0, use_non_chord_patts: bool = False) -> list[GriddedChord]:
        """Returns all patterns from one to up to size points that can be gridded on the tiling
        (only uses active cells)"""
        all_chords = []
        for num_chords in range(1, size//2 + 1):
            all_chords += list(Chord.of_length(num_chords))
        
        chords_on_tiling = set() # this is a set so the non chord diagram pattern code works nicely
        cells = self.active_cells
            
        for chord in all_chords:
            chords_on_tiling.update(filter(self.contains, GriddedChord.all_grids(chord, cells)))

        if use_non_chord_patts:
            bigger_chords = []
            for num_chords in range(size//2 + 1, size + 1):
                bigger_chords += list(map(lambda chord: chord.get_pattern(), Chord.of_length(num_chords)))
                
            for chord_size in range(1, size + 1):
                for chord in bigger_chords:
                    for subchord in list(combinations(chord, chord_size)):
                        subchord = Chord.to_standard(subchord, True)
                        chords_on_tiling.update(filter(self.contains, GriddedChord.all_grids(subchord, cells)))
                    
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
        
        # if we have a req_list of length 0, we are requiring that we have nothing, so we have the empty set
        if any(len(req_list) == 0 for req_list in self.requirements) and len(self.requirements) >= 1:
            return True
        
        sum_max_reqs = 0
        for req_list in self.requirements:
            #if len(req_list) == 0:
                #return True
            req_lengths = [len(req) for req in req_list]
            sum_max_reqs += max(req_lengths)

        #if sum_max_reqs == 0: # then there are no requirements, but we still need a chord of size one
            #sum_max_reqs += 1
        
        # proved maximum size of smallest chord that can be gridded:
        max_len = sum_max_reqs * 2 - 1

        if max_len == -1:
            max_len = 0
       
        all_chords = []
        for num_chords in range(0, max_len + 1):
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

    def add_obstructions(self, gcs: Iterable[GriddedChord], simplify: bool = True, expand: bool = True) -> "Tiling":
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

        active_cells_str = "\nActive cells: "
        for cell in self.active_cells:
            active_cells_str += str(cell)
            active_cells_str += ", "

        return(dimensions + obs_str + reqs_str + link_str + active_cells_str)
    
    def pretty_print(self) -> None:
        """Pretty print the tiling with distinct labels for obstructions and
        requirements. Uses colored terminal output for quick view and also
        provides a LaTeX/TikZ exporter via `pretty_print_latex`.
        """
        # ANSI color codes
        RED = '\033[91m'
        BLUE = '\033[94m'
        YELLOW = '\033[93m'
        RESET = '\033[0m'
        BOLD = '\033[1m'

        dim_x, dim_y = self._cached_properties['dimensions']

        # Assign unique labels to obstructions (A, B, C, ...) and requirements (a, b, c...)
        obs_labels: Dict[GriddedChord, str] = {}
        req_labels: Dict[GriddedChord, str] = {}
        for i, ob in enumerate(self._obstructions):
            obs_labels[ob] = chr(ord('A') + (i % 26)) + (str(i//26) if i//26 else '')
        all_reqs = [r for req_list in self._requirements for r in req_list]
        for i, rq in enumerate(all_reqs):
            req_labels[rq] = chr(ord('a') + (i % 26)) + (str(i//26) if i//26 else '')

        # Build quick lookup: for each cell pick the label of the first obstruction
        cell_marker: Dict[Cell, Tuple[str, str]] = {}  # cell -> (label, color)
        for ob, lab in obs_labels.items():
            for cell in ob._cells:
                cell_marker[cell] = (lab, RED)
        for rq, lab in req_labels.items():
            for cell in rq._cells:
                # only place requirement marker if no obstruction label present
                if cell not in cell_marker:
                    cell_marker[cell] = (lab, BLUE)

        # Header
        print(f"{BOLD}Dimensions: ({dim_x}, {dim_y}){RESET}\n")

        # Draw grid (top to bottom)
        for y in range(dim_y - 1, -1, -1):
            # horizontal border
            print("+", end="")
            for x in range(dim_x):
                print("-+", end="")
            print()
            # contents
            print("|", end="")
            for x in range(dim_x):
                cell = (x, y)
                if cell in cell_marker:
                    lab, color = cell_marker[cell]
                    # short label (max 3 chars) to keep grid tidy
                    short = lab[:3]
                    if color == RED:
                        print(f"{RED}{short}{RESET}|", end="")
                    else:
                        print(f"{BLUE}{short}{RESET}|", end="")
                else:
                    print("  |", end="")
            print()
        # bottom border
        print("+", end="")
        for x in range(dim_x):
            print("-+", end="")
        print("\n")

        # Legend for obstructions
        if obs_labels:
            print(f"{BOLD}{YELLOW}Obstruction legend (red):{RESET}")
            for ob, lab in obs_labels.items():
                print(f"  {RED}{lab}{RESET}: {str(ob)}")
        else:
            print(f"{BOLD}Obstructions:{RESET} (none)")
        print()

        # Legend for requirements
        if req_labels:
            print(f"{BOLD}{YELLOW}Requirement legend (blue):{RESET}")
            for rq, lab in req_labels.items():
                print(f"  {BLUE}{lab}{RESET}: {str(rq)}")
        else:
            print(f"{BOLD}Requirements:{RESET} (none)")
        print()

        # show linkages and assumptions succinctly
        if self._linkages:
            print(f"{BOLD}Linkages:{RESET}")
            for linkage in self._linkages:
                cells = ", ".join(str(coord) for coord in linkage)
                print(f"  {{{cells}}}")
            print()
        if self._assumptions:
            print(f"{BOLD}Assumptions:{RESET}")
            for i, assumption in enumerate(self._assumptions):
                print(f"  {i+1}. {str(assumption)}")
            print()

    def pretty_print_latex(self, filename: str = "tiling_visual.tex", compile_pdf: bool = True) -> str:
        """Export a TikZ picture of the tiling to `filename`.
        If `compile_pdf` is True, attempt to run `pdflatex` (from PATH) to
        produce a PDF next to the .tex file. Returns path to .tex file.
        """
        dim_x, dim_y = self._cached_properties['dimensions']

        # label assignment consistent with pretty_print
        obs_labels = {ob: chr(ord('A') + i % 26) + (str(i//26) if i//26 else '') for i, ob in enumerate(self._obstructions)}
        all_reqs = [r for req_list in self._requirements for r in req_list]
        req_labels = {rq: chr(ord('a') + i % 26) + (str(i//26) if i//26 else '') for i, rq in enumerate(all_reqs)}

        # Calculate scale factor based on number of obstructions and requirements
        total_constraints = len(self._obstructions) + len(all_reqs)
        scale_factor = 1.0 + max(0, (total_constraints - 3) * 0.15)

        tex_lines = []
        tex_lines.append("\\documentclass{standalone}")
        tex_lines.append("\\usepackage{tikz}")
        tex_lines.append("\\begin{document}")
        tex_lines.append(f"\\begin{{tikzpicture}}[scale={scale_factor}]")

        # draw rounded cell boxes
        tex_lines.append("  % draw rounded cell boxes")
        for x in range(dim_x):
            for y in range(dim_y):
                tex_lines.append(f"  \\draw[rounded corners=3pt] ({x},{y}) rectangle ({x+1},{y+1});")

        # helper: cell-local coordinate
        def cell_coord(cell, ux, uy):
            return (cell[0] + ux, cell[1] + uy)

        # draw obstructions: place endpoints inside each cell and draw connections
        for idx, ob in enumerate(self._obstructions):
            col = "red"
            lab = obs_labels[ob]
            
            # Determine vertical levels for each unique chord ID (bottom to top: 0 below 1 below 2, etc)
            unique_chords = sorted(set(ob._chord_dict.keys()))
            chord_to_level = {}
            
            # Distribute obstruction base levels across the vertical space
            num_obs = len(self._obstructions)
            if num_obs == 1:
                obs_base_uy = 0.5
            else:
                # Spread obstructions evenly from 0.1 to 0.9
                obs_base_uy = 0.1 + 0.8 * (idx / max(1, num_obs - 1))
            
            if len(unique_chords) == 1:
                # Single chord: use the obstruction's base level
                chord_to_level[unique_chords[0]] = obs_base_uy
            else:
                # Multiple chords: distribute around the obstruction's base level
                num_chords = len(unique_chords)
                for i, chord_id in enumerate(unique_chords):
                    # Spread chords evenly: chord 0 at lowest, chord N-1 at highest
                    uy = obs_base_uy - 0.1 + (i / max(1, num_chords - 1)) * 0.2 if num_chords > 1 else obs_base_uy
                    chord_to_level[chord_id] = max(0.05, min(0.95, uy))
            
            # Place endpoints left-to-right in the order they appear in the pattern
            # Each endpoint's horizontal position is based on its position in the overall pattern
            pts: Dict[int, Tuple[float, float]] = {}
            num_endpoints = len(ob._pos)
            
            for endpoint_idx in range(num_endpoints):
                cell = ob._pos[endpoint_idx]
                
                # Find which chord uses this endpoint to determine vertical level
                uy = 0.5  # default
                for chord_id, (i1, i2) in ob._chord_dict.items():
                    if endpoint_idx == i1 or endpoint_idx == i2:
                        uy = chord_to_level[chord_id]
                        break
                
                # Horizontal position: left-to-right based on endpoint index
                ux = 0.2 + 0.6 * (endpoint_idx / max(1, num_endpoints - 1)) if num_endpoints > 1 else 0.5
                
                pts[endpoint_idx] = cell_coord(cell, ux, uy)
            
            # draw circles for all endpoints
            for i, (x, y) in pts.items():
                tex_lines.append(f"  \\fill[{col}] ({x},{y}) circle (0.06);")

            # draw polyline connecting all endpoints
            visited = set()
            polyline_points = []
            for idx, chord_id in enumerate(ob._patt):
                if chord_id not in visited:
                    visited.add(chord_id)
                    i1, i2 = ob._chord_dict[chord_id]
                    polyline_points.append(i1)
                    polyline_points.append(i2)
            if len(polyline_points) >= 2:
                path = " -- ".join(f"({pts[i][0]},{pts[i][1]})" for i in polyline_points)
                tex_lines.append(f"  \\draw[{col}, line width=1.2pt] {path};")

        # draw requirements similarly but in blue dashed style
        for idx, rq in enumerate(all_reqs):
            col = "blue"
            lab = req_labels[rq]
            
            # Determine vertical levels for each unique chord ID (bottom to top: 0 below 1 below 2, etc)
            unique_chords = sorted(set(rq._chord_dict.keys()))
            chord_to_level = {}
            
            # Distribute requirement base levels across the vertical space
            num_reqs = len(all_reqs)
            if num_reqs == 1:
                req_base_uy = 0.5
            else:
                # Spread requirements evenly from 0.1 to 0.9
                req_base_uy = 0.1 + 0.8 * (idx / max(1, num_reqs - 1))
            
            if len(unique_chords) == 1:
                # Single chord: use the requirement's base level
                chord_to_level[unique_chords[0]] = req_base_uy
            else:
                # Multiple chords: distribute around the requirement's base level
                num_chords = len(unique_chords)
                for i, chord_id in enumerate(unique_chords):
                    # Spread chords evenly: chord 0 at lowest, chord N-1 at highest
                    uy = req_base_uy - 0.1 + (i / max(1, num_chords - 1)) * 0.2 if num_chords > 1 else req_base_uy
                    chord_to_level[chord_id] = max(0.05, min(0.95, uy))
            
            # Place endpoints left-to-right in the order they appear in the pattern
            # Each endpoint's horizontal position is based on its position in the overall pattern
            pts: Dict[int, Tuple[float, float]] = {}
            num_endpoints = len(rq._pos)
            
            for endpoint_idx in range(num_endpoints):
                cell = rq._pos[endpoint_idx]
                
                # Find which chord uses this endpoint to determine vertical level
                uy = 0.5  # default
                for chord_id, (i1, i2) in rq._chord_dict.items():
                    if endpoint_idx == i1 or endpoint_idx == i2:
                        uy = chord_to_level[chord_id]
                        break
                
                # Horizontal position: left-to-right based on endpoint index
                ux = 0.2 + 0.6 * (endpoint_idx / max(1, num_endpoints - 1)) if num_endpoints > 1 else 0.5
                
                pts[endpoint_idx] = cell_coord(cell, ux, uy)
            
            # draw circles for all endpoints
            for i, (x, y) in pts.items():
                tex_lines.append(f"  \\fill[{col}] ({x},{y}) circle (0.05);")

            # draw polyline connecting all endpoints
            visited = set()
            polyline_points = []
            for idx_pat, chord_id in enumerate(rq._patt):
                if chord_id not in visited:
                    visited.add(chord_id)
                    i1, i2 = rq._chord_dict[chord_id]
                    polyline_points.append(i1)
                    polyline_points.append(i2)
            if len(polyline_points) >= 2:
                path = " -- ".join(f"({pts[i][0]},{pts[i][1]})" for i in polyline_points)
                tex_lines.append(f"  \\draw[{col}, dashed] {path};")

        tex_lines.append("\\end{tikzpicture}")
        tex_lines.append("\\end{document}")

        with open(filename, "w") as f:
            f.write("\n".join(tex_lines))

        # attempt to compile
        if compile_pdf:
            import shutil, subprocess, os
            pdflatex = shutil.which("pdflatex")
            if pdflatex:
                try:
                    subprocess.run([pdflatex, "-interaction=nonstopmode", filename], cwd=os.path.dirname(os.path.abspath(filename)) or ".", check=True, stdout=subprocess.DEVNULL)
                    # open the PDF on macOS
                    pdf_file = filename.replace(".tex", ".pdf")
                    if os.path.exists(pdf_file):
                        subprocess.run(["open", pdf_file], check=False)
                except Exception:
                    pass

        return filename
    
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
    
gc_single_00_00 = GriddedChord(Chord((0, 0)), ((0, 0), (0, 0)))
t_no_restrictions = Tiling((), (), (), (), derive_empty=False)
assert t_no_restrictions.contains(gc_single_00_00)


tiling2 = Tiling(
        obstructions=[
            GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), ) * 4),
            GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
            GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))
        ],
        requirements=[
        ],
    )
atom = Tiling(obstructions=(GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),
                            GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (0, 0), (0, 0))),
                            GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))),
              requirements=((GriddedChord(Chord((0, 0)), ((0,0), (0,0))),),))

atom.pretty_print_latex("atom.tex")

non_crossing = Tiling(obstructions=(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),))

non_crossing.pretty_print_latex("non_crossing.tex")

theorem303 = Tiling(obstructions=(GriddedChord(Chord((0, 1, 2, 0, 1, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
                                  GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
                                  GriddedChord(Chord((0, 1, 2, 1, 0, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),))

theorem303.pretty_print_latex("theorem303.tex")