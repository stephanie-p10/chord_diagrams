"""Tiling: grid classes defined by chord obstructions and requirements.

`Tiling` is the main combinatorial class in this repository. Informally, a
tiling describes a class of gridded chord diagrams constrained by:

- **obstructions**: gridded chord patterns that must be avoided
- **requirements**: tuples/lists of gridded chord patterns, where each list
  represents an “OR” condition (at least one must occur)
- **linkages**: sets of cells that are treated as linked (used by some
  algorithms/strategies)
- **assumptions**: tracking assumptions used by the specification searcher

The implementation maintains cached derived data such as active/empty cells and
row/column maps; initialization can optionally expand, simplify, and/or remove
empty rows/columns.
"""
import sys
from pathlib import Path

_src_root = Path(__file__).resolve().parents[2]  # src
if str(_src_root) not in sys.path:
    sys.path.insert(0, str(_src_root))

import json
from itertools import chain, filterfalse, product, combinations, combinations_with_replacement, groupby
from typing import (Any, Callable, Dict, FrozenSet, Iterable, Iterator, List, Optional, Set, Tuple)
from collections import Counter, defaultdict

from typing_extensions import TypedDict

from comb_spec_searcher import CombinatorialClass, VerificationStrategy
from comb_spec_searcher.exception import StrategyDoesNotApply
from comb_spec_searcher.typing import Parameters

from src.common.chords import GriddedChord, Chord
from src.common.assumptions import TrackingAssumption
from src.common.latex_exporter import export_tiling_to_latex
from src.algorithms.map import RowColMap
from src.algorithms.simplify import SimplifyObstructionsAndRequirements
from src.algorithms.expansion import Expansion
from tilings.misc import intersection_reduce, union_reduce
from tilings.assumptions import (
    ComponentAssumption,
    SkewComponentAssumption,
    SumComponentAssumption,
)

import time

__all__ = ["Tiling"]

Cell = Tuple[int, int]
ReqList = Tuple[GriddedChord, ...]
CellBasis = Dict[Cell, Tuple[List[Chord], List[Chord]]]
CellFrozenSet = FrozenSet[Cell]
Dimension = Tuple[int, int]
GCTuple = Tuple[GriddedChord, ...]


class Tiling(CombinatorialClass):
    """A grid class of gridded chord diagrams.

    **Coordinate system**: cells are 0-indexed `(x, y)` from the bottom-left,
    where `x` is the column index and `y` the row index.

    **Stored inputs** (normalized to tuples during `__init__`):\n
    - `obstructions`: tuple of `GriddedChord`\n
    - `requirements`: tuple of requirement lists (each requirement list is a tuple of `GriddedChord`)\n
    - `linkages`: tuple of tuples of cells\n
    - `assumptions`: tuple of tracking assumptions\n

    **Cached derived properties** include dimensions, empty/active cells, and
    forward/backward row/column maps used by strategies and algorithms.
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

        # Initialization can be relatively expensive (expansion/simplification),
        # so we track timing in case callers want to profile.
        start_time = time.time()

        super().__init__()
        self._linkages = tuple(tuple(link) for link in linkages)
        self._obstructions = tuple(obstructions)
        self._requirements = tuple(tuple(req_list) for req_list in requirements)
        self._assumptions = tuple(assumptions)
        self._cached_properties = {}
        self._derive_empty = derive_empty

        if "dimensions" not in self._cached_properties:
            self._compute_dimensions()

        # currently defaults cells with no requirements or obsturctions to be empty
        # note: should this be defaulted to allowing other chords? - NO
        if "empty_cells" not in self._cached_properties:
            self._prepare_properties(derive_empty)

        if expand:
            self._expand()

        # for simplification of computation, point obstructions to mark empty cells are added after expand.
        # this should not affect the switch to only describe a tiling with chord diagrams, since the only time
        #   the new point obstructions are used is when dealing with parts of the tiling that are empty. Most 
        #   computations are done only on the active cells.
        self._add_point_obs(self._cached_properties["empty_cells"])
            
        # Simplification removes redundant constraints and cleans up patterns
        # that cannot interact with active cells.
        if simplify:
            self._simplify()

        if remove_empty_rows_and_cols:
            self._remove_empty_rows_and_cols()
 
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
        # Fast method of calculating active cells, assumes the user did not do anything "silly"
        # adds all cells that an obstuction larger than a single point uses
        cells_used = union_reduce(
            set(ob.pos) for ob in self.obstructions if len(ob.patt) > 1
        )
        # adds all cells that are used in a requirement to the active cells.
        cells_used.update(
            *(union_reduce(set(comp.pos) for comp in req) for req in self.requirements)
        )

        # calculates dimensions based on what the max cells used in an obstuction or requirement
        max_row = 0
        max_col = 0
        for cell in cells_used:
            max_col = max(max_col, cell[0])
            max_row = max(max_row, cell[1])
        dimensions = (max_col + 1, max_row + 1)

        # finds the cells that have point obstructions - these should be empty
        point_ob_cells = []
        for ob in self.obstructions:
            if len(ob.patt) == 1:
                point_ob_cells.append(ob.pos[0])

        if derive_empty:
            # if we are assuming cells with no information are empty, the active cells
            # are the cells that get mentioned in an ob or req, minus any point cells
            active_cells = cells_used.difference(point_ob_cells)
            empty_cells = tuple(
                cell
                for cell in product(range(dimensions[0]), range(dimensions[1]))
                if cell not in active_cells
            ) 
        else:
            # if we are assuming all cells are active unless explicitly stated, the only empty
            # cells are the cells that have point obstructions, and the active cells are the complement
            empty_cells = tuple([])
            active_cells = tuple(cell 
                                 for cell in product(range(dimensions[0]), range(dimensions[1]))
                                 if cell not in empty_cells) 

        self._cached_properties["active_cells"] = frozenset(active_cells)
        self._cached_properties["empty_cells"] = frozenset(empty_cells)
        self._cached_properties["dimensions"] = dimensions
    
    def _add_point_obs(self, cells: Tuple[Cell]):
        new_obs = list(self._obstructions)
        for cell in cells:
            new_obs.append(GriddedChord.single_cell(Chord((0,)), cell))

        self._obstructions = tuple(new_obs)
  
    # -------------------------------------------------------------
    # Properties and getters
    # -------------------------------------------------------------
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
        """Row/column map from this tiling to a normalized coordinate system.

        Many operations (notably `remove_empty_rows_and_cols`) transform a tiling
        by collapsing unused rows/columns. The resulting coordinate change is
        represented as a `RowColMap` that can be applied to cells, gridded chords,
        and tracking assumptions.

        This property is computed lazily; requesting it may trigger removal of
        empty rows/columns to populate the cache.
        """
        try:
            return self._cached_properties["forward_map"]
        except KeyError:
            self._remove_empty_rows_and_cols()
            return self._cached_properties["forward_map"]

    @property
    def backward_map(self) -> RowColMap:
        """Inverse of `forward_map` (when the forward map is reversible)."""
        try:
            return self._cached_properties["backward_map"]
        except KeyError:
            backward_map = self.forward_map.reverse()
            self._cached_properties["backward_map"] = backward_map
            return backward_map
    
    @property
    def active_cells(self) -> CellFrozenSet:
        """
        Return the active (non-empty) cells of the tiling.

        A cell is considered *empty* if it is marked by a point obstruction used
        internally to represent derived emptiness. Active cells are those not
        marked empty.
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
    
    @property
    def up_to_point_cells(self) -> CellFrozenSet:
        """ 
        Finds the cells that must have 0 or 1 points in them.
        """
        try:
            return self._cached_properties["up_point_cells"]
        except KeyError:
            # generates all chords on the tiling (we only need to check size 2 subchords, but generating this is tricky) TODO
            chords_to_check = self.chords_of_length(self.maximum_length_of_minimum_gridded_chord() + 1)

            # if a cell has more than two points in it, it must contain the pattern 00, 01, or 10
            # this means that the cells that can have one point are the cells that do not contain 00, 01 or 10
            # in other words, they do not contain 00, or a size two chord pattern with at least two points in that cell
            bad_cells = []

            # finding which cells have 00
            for cell in self.active_cells:
                if self.contains(GriddedChord(Chord((0, 0)), (cell, cell))):
                    bad_cells.append(cell)

            # finding the cells that contain two points of a size two chord
            for gc in chords_to_check:
                pos_list_counter = Counter(gc.pos)
                for cell in pos_list_counter.keys():
                    if pos_list_counter[cell] > 1:
                        bad_cells.append(cell)

            result = (cell for cell in self.active_cells if cell not in bad_cells)
            self._cached_properties["up_to_point_cells"] = result
            return result

    @property
    def point_cells(self) -> CellFrozenSet:
        """ 
        Finds cells that must have exactly one point in them.
        """
        try:
            return self._cached_properties["point_cells"]
        except KeyError:

            point_cells = []
            for cell in self.up_to_point_cells:
                if cell in self.positive_cells:
                    point_cells.append(cell)

            self._cached_properties["point_cells"] = point_cells
            return point_cells
        
    @property
    def chord_cells(self) -> CellFrozenSet:
        """Not tested, but should return cells that have exactly one chord in them"""
        try:
            return self._cached_properties["chord_cells"]
        except KeyError:
            # chord cells are cells that contain exactly one chord (both ends)
            # this means they contain 00, and avoid 01, 10
            possible_00_chords = [GriddedChord(Chord((0, 0)), (cell, cell)) for cell in self.active_cells]
            cells_with_00 = []
            
            # check which cells must contain a 00 chord diagram pattern
            for reqlist in self.requirements:
                if len(reqlist) >= 1:
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
        
    @property
    def point_col_cells(self):
        """Returns point cells that are isolated in their column"""
        try:
            return self._cached_properties["point_col_cells"]
        except KeyError:
            cells = [cell for cell in self.point_cells if self.only_cell_in_col(cell)]
            self._cached_properties["point_col_cells"] = cells
            return cells
        
    @property
    def chord_row_cells(self) -> list[Cell]:
        """Returns cells that contain a chord that is isolated in its own row"""
        try:
            return self._cached_properties["chord_row_cells"]
        except KeyError:
            # finds cells that contain exactly one chord and are the only active cell in their row
            row_isolated_chord_cells = [cell for cell in self.chord_cells if self.only_cell_in_row(cell)]

            # finds cells that contain exatly one point in rows that have exactly one chord
            active_cells_in_row_counter = Counter(cell[1] for cell in self.active_cells)
            point_cells_in_row_counter = Counter(cell[1] for cell in self.point_cells)
            # finds cells that are in a row with exactly point cells, and no additional active cells
            row_isolated_point_cells = [cell 
                                        for cell in self.point_cells 
                                        if point_cells_in_row_counter[cell[1]] == 2 and 
                                        active_cells_in_row_counter[cell[1]] == 2]
            
            all_row_isolated_cells = row_isolated_chord_cells + row_isolated_point_cells
            self._cached_properties["chord_row_cells"] = all_row_isolated_cells
            return all_row_isolated_cells
  
    def sub_tiling(
        self,
        cells: Iterable[Cell],
        factors: bool = False,
        add_assumptions: Iterable[TrackingAssumption] = tuple(),
    ) -> "Tiling":
        """Restrict this tiling to a specified set of cells.

        Only obstructions/requirements/linkages fully supported on `cells` are
        kept (and assumptions are filtered accordingly). This is the primary
        operation used by factorization strategies to produce child tilings.

        - **cells**: iterable of cells to keep.
        - **factors**: optimization flag used by factorization code paths. When
          true, the implementation may use the first position as a quick filter
          before checking all positions.
        - **add_assumptions**: extra assumptions to include in the result.
        """
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

    def maximum_length_of_minimum_gridded_chord(self) -> int:
        """Returns the maximum length of the minimum gridded chord diagram that
        can be gridded on the tiling.
        """
        if not self.requirements:
            return 1 # should this be 0? sTODO
        
        sum_largest_req_lengths = 0
        for reqs in self.requirements:
            if len(reqs) != 0:
                req_lengths = [max(req.patt) + 1 for req in reqs if len(req.patt) > 0]
                sum_largest_req_lengths += max(req_lengths)
        
        return 2 * sum_largest_req_lengths - 1

    def _minimum_size_of_object(self):
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

        if max_len == -1: # then there are no requirements, but we still need a chord of size one
            max_len = 0

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
        
        return -1 # this maybe should be 0?
    
    def minimum_size_of_object(self) -> int:
        try:
            min_size = self._cached_properties["minimum_size_of_object"]
            return min_size
        except KeyError:
            min_size = self._minimum_size_of_object()
            self._cached_properties["minimum_size_of_object"] = min_size
            return min_size
    
    # s TODO fix all_chords_on_tiling to use this
    def chords_of_length(self, length: int, use_non_chord_patts: bool = False) -> List[GriddedChord]:
        """Returns all chord diagrams on the tiling of length length"""
        all_chords = list(Chord.of_length(length))
        cells = list(self.active_cells)
        chords_on_tiling = set()
        for chord in all_chords:
            chords_on_tiling.update(filter(self.contains, GriddedChord.all_grids(chord, cells)))

        return list(chords_on_tiling)

    # sToDo this is inefficient and in permutations there is a class that does this efficiently
    # gridded_perms_of_length -> all_chords_on_tiling
    # s TODO check where this is used, see if the empty chord should be included?
    def all_chords_on_tiling(self, size: int = 0, use_non_chord_patts: bool = False) -> List[GriddedChord]:
        """Enumerate all gridded chord patterns up to a given size.

        This is an expensive, brute-force enumerator used by a few inferral and
        emptiness checks. It only grids over the tiling's active cells.

        - **size**: maximum number of points to consider (patterns are generated
          by chord length, then gridded into active cells).
        - **use_non_chord_patts**: if True, additionally considers certain
          standardized subpatterns derived from larger chord patterns.
        """
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
    
    def build_all_chords_on_tiling(self, size: int = 0) -> List[GriddedChord]:
        # Sort obstructions by size
        obstructions_by_size = dict(groupby(sorted(self.obstructions, key=len), key=len))

        single_chords = [] # All allowable single chords on the tiling
        for cell1, cell2 in combinations_with_replacement(self.active_cells, 2):
            if cell1[1] == cell2[1]:
                if cell1[0] <= cell2[0]:
                    chord = GriddedChord(Chord((0, 0)), (cell1, cell2))
                else:
                    chord = GriddedChord(Chord((0, 0)), (cell2, cell1))
                # Only add the chord if it's valid
                if chord.avoids(*obstructions_by_size.get(1, ())):
                    single_chords.append(chord)

        # All gridded chord diagrams we've found that were valid
        result = []

        def build_from_gc(starting_gc: GriddedChord, 
                          potential_chords_to_insert: List[GriddedChord],
                          ) -> List[GriddedChord]:
            """
            Returns all chord diagrams (with length <= size), including starting_gc,
                that can be built by adding chords from potential_chords_to_insert

            Does not modify input
            """
            
            result.append(starting_gc)
            if len(starting_gc) == size:
                # Don't recurse if we've reached the maximum size
                return
            
            chords_inserted = []
            newly_constructed_gcs = []
            min_source_idx_in_patt = 0 if len(starting_gc) == 0 else 1 + starting_gc.chord_dict[len(starting_gc) - 1][0] # one more than the highest chord in starting gc
            min_source_pos = (0, 0) if len(starting_gc) == 0 else starting_gc.pos[min_source_idx_in_patt - 1] # position of highest chord in starting gc

            for gc_to_add in potential_chords_to_insert:
                # If the source cell above the previous highest source...
                if gc_to_add.pos[0][0] >= min_source_pos[0] and gc_to_add.pos[0][1] >= min_source_pos[1]:
                    
                    min_source_idx_in_col, max_source_idx, min_sink_idx, max_sink_idx = starting_gc.get_bounding_indices(gc_to_add.pos[0][0], gc_to_add.pos[1][0])
                    # minimum index for next chord
                    min_source_idx = max(min_source_idx_in_col, min_source_idx_in_patt)
                    # for each source index and sink index after the source, from the minimums to the maximums inclusive...
                    for source_idx in range(min_source_idx, max_source_idx + 1):
                        for sink_idx in range(max(source_idx, min_sink_idx), max_sink_idx + 1):
                            # get the chord diagram where we've added this chord
                            extended_gc = starting_gc.insert_specific_chord(gc_to_add.pos[0][1], 
                                                                         gc_to_add.pos[0][0], 
                                                                         gc_to_add.pos[1][0],
                                                                         source_idx,
                                                                         sink_idx)
                            if extended_gc != None and extended_gc.avoids(*self.obstructions):  # this check is not optimized
                                chords_inserted.append(gc_to_add)
                                newly_constructed_gcs.append(extended_gc)

            for gc_to_build_on in newly_constructed_gcs:
                build_from_gc(gc_to_build_on, chords_inserted)

        build_from_gc(GriddedChord.empty_chord(), single_chords)
        return [gc for gc in result if all(gc.contains(*req) for req in self._requirements)]
      
    def objects_of_size(self, size):
        chords = Chord.of_length(size)

        grids = []
        for chord in chords:
            for grid in GriddedChord.all_grids(chord):
                if self.contains(grid):
                    grids.append(grid)

        return grids

    def is_atom(self):
        """Return True if the Tiling is a single gridded chord."""

        min_size = self.minimum_size_of_object()
        if min_size == -1:
            return False

        grids_of_minimum_size = self.chords_of_length(min_size)
        if len(grids_of_minimum_size) != 1:
            return False
        
        grids_of_bigger_size = self.chords_of_length(min_size + 1)
        if len(grids_of_bigger_size) > 0:
            return False

        return True
    
    def is_empty(self) -> bool:
        """Checks if the tiling is empty.

        Tiling is empty if it has been inferred to be contradictory due to
        contradicting requirements and obstructions or no gridded chord
        can be gridded on the tiling.
        """
        # if any obstruction is empty, no chords can be gridded since all grids contain the empty grid.
        if any(ob.is_empty() for ob in self.obstructions):
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
                if self.contains(gc):
                    return False

        return True

    def contains(self, gc: GriddedChord) -> bool:
        """Return True if `gc` is valid for this tiling.

        A gridded chord diagram is valid if it:
        - satisfies every requirement list (contains at least one requirement in
          each list),
        - avoids every obstruction, and
        - respects all linkage connectivity constraints.
        """
        has_reqs = all(gc.contains(*req) for req in self._requirements)
        avoids_ob = not any(gc.contains(ob) for ob in self._obstructions)
        links_connected = all(gc.is_connected(cells) for cells in self._linkages)
        #print(has_reqs, avoids_ob, links_connected)
        return has_reqs and avoids_ob and links_connected
     
    # -------------------------------------------------------------
    # Not used anywhere?
    # -------------------------------------------------------------
    @property
    def required_chords(self) -> List[GriddedChord]:
        """Returns a list of chords that must be contained in any valid gridded chord diagram of the Tiling.
        
        A chord is required if it's contained in all requirements of at least one requirement list.
        Since at least one requirement from each list must be satisfied, any chord appearing 
        in all requirements of a list must appear in any valid gridding.
        """
        try:
            return self._cached_properties["required_chords"]
        except KeyError:
            required_set = set()
            
            for req_list in self.requirements:
                if len(req_list) > 0:

                    # Determine maximum chord size to check based on requirement patterns
                    max_chord_size = 2 * max((max(req.patt) + 1 for req in req_list if len(req.patt) > 0), default=1)
                
                    # Generate candidate chords and check which ones appear in all requirements
                    for length in range(1, max_chord_size):
                        for chord in Chord.of_length(length):
                            # A chord is required if its pattern is contained in ALL requirements' patterns of this list
                            if all(req.patt != () and chord.get_pattern() in req.patt for req in req_list):
                                required_set.add(chord)
            
            required_chords_list = list(required_set)
            self._cached_properties["required_chords"] = required_chords_list
            return required_chords_list
    
        # not used anywhere?

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
  
    # not used anywhere?
    def contained_requirements(self, requirement: GriddedChord) -> List[GriddedChord]:
        """Find all sub-requirements of the given requirement. This maintains the 
           gridding of the requirement.
        """
        # number of chords (labels 0..n-1)
        n = len(requirement._chord)
        if n == 0:
            return []

        result: List[GriddedChord] = []
        seen: set = set()
        labels = range(n)

        for r in range(1, n + 1):
            for subset in combinations(labels, r):
                subset_set = set(subset)
                # build pattern and corresponding positions preserving original order
                patt = [val for val in requirement._patt if val in subset_set]
                positions = [pos for idx, pos in enumerate(requirement._pos) if requirement._patt[idx] in subset_set]
                new_chord = Chord.to_standard(patt)
                key = (tuple(new_chord.get_pattern()), tuple(positions))
                if key in seen:
                    continue
                seen.add(key)
                result.append(GriddedChord(new_chord, positions))

        return result

    # -------------------------------------------------------------
    # Algorithms
    # -------------------------------------------------------------
    def _expand(self) -> None:
        expansion_class = Expansion(self._obstructions, self._requirements, self.dimensions, self.active_cells)
        expansion_class.expand_obstructions()
        expansion_class.expand_requirements()
        self._obstructions = tuple(expansion_class._obstructions)
        self._requirements = tuple(expansion_class._requirements)

    def _simplify(self) -> None:
        simplify_algo = SimplifyObstructionsAndRequirements(self.obstructions, 
                                                            self.requirements, 
                                                            self._cached_properties["dimensions"], 
                                                            self._cached_properties["active_cells"], 
                                                            self._cached_properties["empty_cells"])
        simplify_algo.simplify()
        self._obstructions = simplify_algo.obstructions
        self._requirements = simplify_algo.requirements
        self._cached_properties["active_cells"] = simplify_algo.active_cells
        self._cached_properties["empty_cells"] = simplify_algo.empty_cells

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
            self._obstructions = (GriddedChord.single_cell(Chord((0,0)), (0, 0)),)
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

    # -------------------------------------------------------------
    # Cell methods
    # -------------------------------------------------------------
    def cells_in_row(self, row: int) -> CellFrozenSet:
        """Return all active cells in row."""
        return frozenset((x, y) for (x, y) in self.active_cells if y == row)

    def cells_in_col(self, col: int) -> CellFrozenSet:
        """Return all active cells in column."""
        return frozenset((x, y) for (x, y) in self.active_cells if x == col)
    
    def cell_ob_basis(self) -> Dict[Cell, List[Chord]]:
        """Returns a dictionary from cells to basis.

        The ob basis for each cell is a list of chord diagrams containing the 
        patterns of the obstructions localized in the cell.
        """
        try:
            return self._cached_properties["cell_ob_basis"]
        except KeyError:
            obdict: Dict[Cell, List[Chord]] = defaultdict(list)
            for ob in self.obstructions:
                if ob.is_localized():
                    cell = ob.pos[0]
                    obdict[cell].append(ob.patt)

            # puts point obs in any cells that are inferred to be empty
            for cell in self._cached_properties["empty_cells"]:
                obdict[cell] = [Chord((0,))]

            all_cells = product(range(self.dimensions[0]), range(self.dimensions[1]))
            resdict = {cell: obdict[cell] for cell in all_cells}
            self._cached_properties["cell_ob_basis"] = resdict
            return resdict

    def add_list_requirement(self, req_list: Iterable[GriddedChord]) -> "Tiling":
        """
        Return a new tiling with an additional requirement list.

        Each requirement list represents an OR-condition: an object in the
        tiling must contain at least one of the patterns in the list.
        """
        new_req = tuple(sorted(req_list))
        return Tiling(
            self._obstructions,
            sorted(self._requirements + (new_req,)),
            self._linkages,
            self._assumptions,
        )

    def add_requirement(self, patt: Chord, pos: Iterable[Cell]) -> "Tiling":
        """Convenience wrapper to add a single requirement as a 1-item list."""
        new_req_list = (GriddedChord(patt, pos),)
        return self.add_list_requirement(new_req_list)

    def add_obstructions(self, gcs: Iterable[GriddedChord], simplify: bool = True, expand: bool = True) -> "Tiling":
        """Return a new tiling with additional obstructions.

        - **simplify**: whether to run simplification after adding.
        - **expand**: whether to run expansion after adding.
        """
        new_obs = tuple(gcs)
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
    
    def obs_in_cell(self, cell: Cell) -> Tuple[int, ...]:
        """Returns the indices obstructions in cell"""

    def reqs_in_cell(self, cell: Cell) -> Tuple[int, ...]:
        """Returns the indices of requirments lists where at least one requirment in the list has a point in cell"""

    # -------------------------------------------------------------
    # Dunder methods
    # -------------------------------------------------------------
    def __hash__(self) -> int:
        return (
            hash(self._requirements)
            ^ hash(self._obstructions)
            ^ hash(self._assumptions)
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Tiling):
            return False
        #print(self.obstructions == other.obstructions)
        #print(self.requirements == other.requirements)
        #print(self.linkages == other.linkages)
        #print(self.assumptions == other.assumptions)
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

    def __contains__(self, gp: GriddedChord) -> bool:
        """Test if a gridded chord is griddable on the given tiling."""
        return (
            gp.avoids(*self.obstructions)
            and all(gp.contains(*req) for req in self.requirements)
            and all(
                (len(linkage) == 0) or gp.is_connected(list(linkage))
                for linkage in self.linkages
            )
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

    def __str__(self) -> str:
        
        # pylint: disable=too-many-locals
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
        for cell, obstructions in sorted(self.cell_ob_basis().items()):
            basis = list(sorted(obstructions))
            if basis == [Chord((0,))]:
                continue # don't print anything if the cell is empty
            # the block is the basis and whether or not positive
            block = (tuple(basis), cell in self.positive_cells)
            label = labels.get(block)
            if label is None:
                if cell in self.up_to_point_cells:
                    if cell in self.positive_cells:
                        label = "\u25cf"
                        block = ((Chord((0, 1)), Chord((1, 0)), Chord((0, 0))), True)
                    else:
                        label = "\u25cb"
                        block = ((Chord((0, 1)), Chord((1, 0)), Chord((0, 0))), False)
                else:
                    label = str(curr_label)
                    curr_label += 1
                labels[block] = label
            row_index_from_top = dim_j - cell[1] - 1
            index = (2 * row_index_from_top + 1) * row_width + 2 * cell[0] + 1
            result[index] = label

        '''chord_rows_dict = {}
        for col, row in self.chord_row_cells:
            try:
                first_col = chord_rows_dict[row][0]
                if first_col <= col:
                    chord_rows_dict[row] = [first_col, col]
                else:
                    chord_rows_dict[row] = [col, first_col]
            except KeyError:
                chord_rows_dict[row] = [col]

        for row, cols in chord_rows_dict.items():
            for col in range(cols[0]+1, cols[1]):
                row_index_from_top = dim_j - row - 1
                index = (2 * row_index_from_top + 1) * row_width + 2 * col + 1
                result[index] = "\u2500"''' # prints lines to represent chords on any cells that are chord cells

        # Legend at bottom
        for block, label in sorted(labels.items(), key=lambda x: x[1]):
            basis_el, positive = block
            result.append(label)
            result.append(": ")
            if basis_el == (Chord((0, 1)), Chord((1, 0)), Chord((0, 0))) and positive:
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

        return "".join(result)
    
        '''OLD CODE: does not print grid, just lists attributes
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

        return(dimensions + obs_str + reqs_str + link_str + active_cells_str)'''

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
    
    # -------------------------------------------------------------
    # Pretty printing
    # -------------------------------------------------------------
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
        return export_tiling_to_latex(self, filename, compile_pdf)
    
    def to_html_representation(self) -> str:
        """Returns an html representation of the tilings object"""
        # pylint: disable=too-many-locals
        # stylesheet for tiling
        style = """
            border: 1px solid;
            width: 24px;
            height: 24px;
            text-align: center;
            """
        dim_i, dim_j = self.dimensions
        result = []
        # Create tiling html table
        result.append("<table> ")
        for _ in range(dim_j):
            result.append("<tr>")
            for _ in range(dim_i):
                result.append(f"<th style='{style}'>")
                result.append(" ")
                result.append("</th>")
            result.append("</tr>")
        result.append("</table>")
        labels: Dict[Tuple[Tuple[Chord, ...], bool], str] = {}

        # Put the sets in the tiles

        # How many characters are in a row in the grid
        row_width = 3 * dim_i + 2
        curr_label = 1
        for cell, obstructions in sorted(self.cell_ob_basis().items()):
            basis = list(sorted(obstructions))
            if basis == [Chord((0,))]:
                continue
            # the block, is the basis and whether or not positive
            block = (tuple(basis), cell in self.positive_cells)
            label = labels.get(block)
            if label is None:
                if cell in self.up_to_point_cells:
                    if cell in self.positive_cells:
                        label = "\u25cf"
                        block = ((Chord((0, 1)), Chord((1, 0)), Chord((0, 0))), True)
                    else:
                        label = "\u25cb"
                        block = ((Chord((0, 1)), Chord((1, 0)), Chord((0, 0))), False)
                else:
                    label = str(curr_label)
                    curr_label += 1
                labels[block] = label
            row_index_from_top = dim_j - cell[1] - 1
            index = row_index_from_top * row_width + cell[0] * 3 + 3
            result[index] = label

        '''chord_rows_dict = {}
        for col, row in self.chord_row_cells:
            try:
                first_col = chord_rows_dict[row][0]
                if first_col <= col:
                    chord_rows_dict[row] = [first_col, col]
                else:
                    chord_rows_dict[row] = [col, first_col]
            except KeyError:
                chord_rows_dict[row] = [col]

        for row, cols in chord_rows_dict.items():
            for col in range(cols[0]+1, cols[1]):
                row_index_from_top = dim_j - row - 1
                index = row_index_from_top * row_width + col * 3 + 3
                result[index] = "\u2500"''' #puts lines for chords between point cells

        # adds background color in cells where assumption happens
        #result = self._handle_html_assumption(result, style)
        return "".join(result)



if __name__ == "__main__":
    gc_single_00_00 = GriddedChord(Chord((0, 0)), ((0, 0), (0, 0)))
    t_no_restrictions = Tiling((), (), (), (), derive_empty=False)
    assert t_no_restrictions.contains(gc_single_00_00)

