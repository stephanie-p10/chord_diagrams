import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from itertools import chain, filterfalse, product
from typing import TYPE_CHECKING, Dict, FrozenSet, Iterable, List, Optional, Tuple
from collections import Counter

from misc import DIR_EAST, DIR_NONE, DIR_NORTH, DIR_SOUTH, DIR_WEST, DIRS
from assumptions import TrackingAssumption
from chords import Chord, GriddedChord
from tiling import Tiling

Cell = Tuple[int, int]
Dir = int
ListRequirement = List[GriddedChord]
Linkages = List[Cell]
LinkCache = Dict[Cell, List[Linkages]]
ObsCache = Dict[Cell, List[GriddedChord]]
ReqsCache = Dict[Cell, List[ListRequirement]]
AssumpCache = Dict[Cell, List[TrackingAssumption]]

'''
    # used for translating cells in linkages
    def cell_translation(self, cell_to_move: Cell, cell_placed: Cell) -> Cell:
        """
        Return how cell_to_move is stretched assuming that a point is placed into cell_placed.
        """
        x, y = cell_to_move
        return (
            cell_to_move[0] + 2 if self.own_col and cell_to_move[0] >= cell_placed[0] else x,
            cell_to_move[1] + 2 if (self.own_row and cell_to_move[1] >= cell_placed[1]) else y,
        )

    def stretched_obstructions(self, cell: Cell) -> List[GriddedChord]:
        """
        Return all of the stretched obstructions that are created if placing a
        point in the given cell.
        """
        if cell not in self._stretched_obstructions_cache:
            self._stretched_obstructions_cache[cell] = self._stretch_gridded_chords(
                self._tiling.obstructions, cell
            )
        return self._stretched_obstructions_cache[cell]

    def stretched_requirements(self, cell: Cell) -> List[ListRequirement]:
        """
        Return all of the stretched requirements that are created if placing a
        point in the given cell.
        """
        if cell not in self._stretched_requirements_cache:
            self._stretched_requirements_cache[cell] = [
                self._stretch_gridded_chords(req_list, cell)
                for req_list in self._tiling.requirements
            ]
        return self._stretched_requirements_cache[cell]
    
    def stretched_linkages(self, cell: Cell) -> List[Linkages]:
        """
        Retrun the stretched linkages obtained from placing a point in cell
        """
        if cell not in self._stretched_linkages_cache:
            self._stretched_linkages_cache[cell] = [
                [self.cell_translation(c) for c in link] 
                for link in self._tiling.linkages]
        return self._stretched_linkages_cache[cell]

    def stretched_assumptions(self, cell: Cell) -> List[TrackingAssumption]:
        """
        Return all of the stretched assumptions that are created if placing a
        point in the given cell.
        """
        if cell not in self._stretched_assumptions_cache:
            self._stretched_assumptions_cache[cell] = [
                assump.__class__(self._stretch_gridded_chords(assump.gcs, cell))
                for assump in self._tiling.assumptions
            ]
        return self._stretched_assumptions_cache[cell]

    # sCN: _stretched_obstructions_requirements_and_assumptions to _stretched_obstructions_requirements_linkages_and_assumptions
    def _stretched_obstructions_requirements_linkages_and_assumptions(
        self, cell: Cell
    ) -> Tuple[List[GriddedChord], List[ListRequirement], List[TrackingAssumption]]:
        """
        Return all of the stretched obstruction and requirements assuming that
        a point is placed in cell.
        """
        stretched_obs = self.stretched_obstructions(cell)
        stretched_reqs = self.stretched_requirements(cell)
        stretched_assump = self.stretched_assumptions(cell)
        stretched_links = self.stretched_linkages(cell)
        point_obs = self._point_obstructions(cell)
        point_req = self._point_requirements(cell)
        return stretched_obs + point_obs, stretched_reqs + point_req, stretched_links, stretched_assump

    @staticmethod
    def _farther(c1: Cell, c2: Cell, direction: Dir) -> bool:
        """Return True if c1 is farther in the given direction than c2."""
        if direction == DIR_EAST:
            return c1[0] > c2[0]
        if direction == DIR_WEST:
            return c1[0] < c2[0]
        if direction == DIR_NORTH:
            return c1[1] > c2[1]
        if direction == DIR_SOUTH:
            return c1[1] < c2[1]
        raise Exception("Invalid direction")

    def forced_obstructions_from_requirement(
        self,
        gcs: Iterable[GriddedChord],
        indices: Iterable[int],
        cell: Cell,
        direction: Dir,
    ) -> List[GriddedChord]:
        """
        Return the obstructions required to ensure that the placed point is
        the direction most occurence of a point used in an occurrence of any
        gridded chord diagrams in gcs.

        In particular, this returns the list of obstructions that are stretched
        from any gridded chord diagram in gcs in which the point at idx is
        farther in the given direction than the placed cell.

        Tells you what you can't have for the placed point to be the directionmost
        """
        placed_cell = self._placed_cell(cell)
        res = []
        if cell in self._tiling.point_cells:
            x, y = placed_cell
            if self.own_row:
                res.append(GriddedChord.point_chord((x, y + 1)))
                res.append(GriddedChord.point_chord((x, y - 1)))
            if self.own_col:
                res.append(GriddedChord.point_chord((x + 1, y)))
                res.append(GriddedChord.point_chord((x - 1, y)))

        for idx, gc in zip(indices, gcs):
            # if cell is farther in the direction than gc[idx], then don't need
            # to avoid any of the stretched grided perms
            if not self._farther(cell, gc.pos[idx], direction):
                for stretched_gc in self._stretch_gridded_chord(gc, cell):
                    if self._farther(stretched_gc.pos[idx], placed_cell, direction):
                        res.append(stretched_gc)
        return res

    def _remaining_requirement_from_requirement(
        self, gcs: Iterable[GriddedChord], indices: Iterable[int], cell: Cell
    ) -> List[GriddedChord]:
        """
        Return the requirements required to ensure that the placed point can be
        extended to be direction most occurrece of a point at the index of the
        gc in gcs.

        In particular, this returns the requirements that come from stretching
        a gridded chord diagram in gcs, such that the point at idx is the placed
        cell.

        What you must have for the point at idx to be placed
        """
        placed_cell = self._placed_cell(cell)
        res = []
        for idx, gc in zip(indices, gcs):
            if gc.pos[idx] == cell:
                for stretched_gp in self._stretch_gridded_chord(gc, cell):
                    if stretched_gp.pos[idx] == placed_cell:
                        res.append(stretched_gp)
        return res

    def place_point_of_gridded_permutation(
        self, gc: GriddedChord, idx: int, direction: Dir
    ) -> "Tiling":
        """
        Return the tiling where the placed point corresponds to the
        directionmost (the furtest in the given direction, ex: leftmost point)
        occurrence of the point at idx in gc.
        """
        return self.place_point_of_req((gc,), (idx,), direction)[0]

    def place_point_of_req(
        self,
        gcs: Iterable[GriddedChord],
        indices: Iterable[int],
        direction: Dir,
        include_not: bool = False,
        cells: Optional[Iterable[Cell]] = None,
    ) -> Tuple["Tiling", ...]:
        """
        Return the tilings where the placed point corresponds to the directionmost
        (the furtest in the given direction, ex: leftmost point) of an occurrence
        of any point [idx, gc(idx)] for gridded chord diagrams in gcs, and idx in indices
        """
        if cells is not None:
            cells = frozenset(cells)
        else:
            cells = frozenset(gc.pos[idx] for idx, gc in zip(indices, gcs))
        res = []
        for cell in sorted(cells):
            stretched = self._stretched_obstructions_requirements_linkages_and_assumptions(cell)
            (obs, reqs, links, assump) = stretched
            remaining_req = self._remaining_requirement_from_requirement(gcs, indices, cell)

            if direction == DIR_NONE:
                res.append(self._tiling.__class__(obs, reqs + [remaining_req], links, assump))
                if include_not:
                    res.append(self._tiling.__class__(obs + remaining_req, reqs, links, assump))
                continue
            forced_obs = self.forced_obstructions_from_requirement(
                gcs, indices, cell, direction
            )
            forced_obs = [
                o1
                for o1 in forced_obs
                # if no other obstruction is contained in o1
                if not any(o2 in o1 for o2 in filterfalse(o1.__eq__, forced_obs))
            ]
            # reduces redundant obstructions
            reduced_obs = [o1 for o1 in obs if not any(o2 in o1 for o2 in forced_obs)]
            # add the forced obstructions as long as they are not redundant
            reduced_obs.extend(filterfalse(reduced_obs.__contains__, forced_obs))
            res.append(
                self._tiling.__class__(
                    reduced_obs,
                    reqs + [remaining_req],
                    links,
                    assumptions=assump
                )
            )
        return tuple(res)

    def place_point_in_cell(self, cell: Cell, direction: Dir) -> "Tiling":
        """
        Return the tiling in which a point is placed in the given direction and
        cell.
        """
        point_req = GriddedChord(Chord(0,), (cell,))
        return self.place_point_of_req((point_req,), (0,), direction)[0]

    def col_placement(self, index: int, direction: Dir) -> Tuple["Tiling", ...]:
        """
        Return the list corresponding the index column being placed in the
        given direction.
        """
        assert direction in (DIR_EAST, DIR_WEST)
        req = [GriddedChord(Chord((0,)), (cell,)) for cell in self._tiling.cells_in_col(index)]
        return self.place_point_of_req(req, tuple(0 for _ in req), direction)

    def row_placement(self, index: int, direction: Dir) -> Tuple["Tiling", ...]:
        """
        Return the list corresponding the index row being placed in the given
        direction.
        """
        assert direction in (DIR_NORTH, DIR_SOUTH)
        req = [GriddedChord(Chord((0,)), (cell,)) for cell in self._tiling.cells_in_row(index)]
        return self.place_point_of_req(req, tuple(0 for _ in req), direction)

    def empty_col(self, index: int) -> "Tiling":
        """
        Return the tiling in which the index row is empty.
        """
        return self._tiling.add_obstructions(
            tuple(
                GriddedChord((0,), (cell,)) for cell in self._tiling.cells_in_col(index)
            )
        )

    def empty_row(self, index: int) -> "Tiling":
        """
        Return the tiling in which the index row is empty.
        """
        return self._tiling.add_obstructions(
            tuple(
                GriddedChord((0,), (cell,)) for cell in self._tiling.cells_in_row(index)
            )
        )
    '''

class RequirementPlacement:
    """
    The requirement placement container class.

    Places points onto own row, own col, or both.

    INPUTS:
        - `tiling`: The tilings to perform the placement with
        - `own_row`: Indiciate to place the point on its own row
        - `own_col`: Indiciate to place the point on its own column
        - `dirs`: The directions used for placement (default to all
          directions).
          The possible directions are:
            - `misc.DIR_NORTH`
            - `misc.DIR_SOUTH`
            - `misc.DIR_EAST`
            - `misc.DIR_WEST`
    """
    #sToDo: update to assume the tiling is described by valid chords. 
    # Maybe check a valid_chords parameter in the tiling that says whether chords are simplified or not?
    def __init__(
        self,
        tiling: "Tiling",
        own_row: bool = True,
        own_col: bool = True,
        dirs: Iterable[int] = tuple(DIRS),
    ):
        if not own_row and not own_col:
            raise ValueError("Must place on own row or on own column.")
        #if own_col and not own_row:
            #raise NotImplementedError
        
        assert all(d in DIRS for d in dirs), "Got an invalid direction"
        self._tiling = tiling
        self.own_row = own_row
        self.own_col = own_col
        if self.own_row and self.own_col:
            self.directions = frozenset(DIRS)
        elif self.own_row:
            self.directions = frozenset((DIR_NORTH, DIR_SOUTH))
        elif self.own_col:
            self.directions = frozenset((DIR_EAST, DIR_WEST))
        self.directions = frozenset(dirs).intersection(self.directions)
        assert self.directions, "No direction to place"

    def chord_already_placed(self, gc: GriddedChord, dir: int) -> bool:
        """
        Determine if the gridded chord diagram gc has already been placed at index idx.
        """
        (idx, val) = self.directionmost_point(gc, dir)
        source_cell = gc.pos[idx]
        sink_cell = gc.pos[gc.chord_dict[val][1]]
        source_placed_in_col = source_cell in self._tiling.point_col_cells
        sink_placed_in_col = sink_cell in self._tiling.point_col_cells
        source_placed_in_row = source_cell in self._tiling.chord_row_cells
        sink_placed_in_row = sink_cell in self._tiling.chord_row_cells

        placed_in_row = source_placed_in_row and sink_placed_in_row
        placed_in_col = source_placed_in_col and sink_placed_in_col

        if self.own_col and self.own_row:
            return placed_in_row and placed_in_col
        if self.own_row:  # Only placing in own row
            return placed_in_row
        if self.own_col:  # Only placing in own column
            return placed_in_col
        raise Exception("Not placing at all!!")
    
    def chords_already_placed(self, gcs: Iterable[GriddedChord], dir: int) -> bool:
        all_placed = all(self.chord_already_placed(gc, dir) for gc in gcs)
        return all_placed
    
    def _point_translation(self, gc: GriddedChord, idx: int, point_placed: Cell) -> tuple[int]:
        """Translates the cell in gc at idx assuming a point is placed at point_placed """
        x, y = gc.pos[idx]
        x = x + 2 if self.own_col and idx >= point_placed[0] else x
        y = y + 2 if (self.own_row and gc.patt[idx] >= point_placed[1]) else y
        return (x, y)

    # this translates a given chord around a placed point that was not in this chord
    # tested!
    def _gridded_chord_translation(
        self, gc: GriddedChord, placed_cell: Cell
    ) -> GriddedChord:
        """
        Return the gridded chord diagram with its positions translated around
        a point placed at placed_cell = (idx, val). Assumes no cell in gc was 
        the cell being placed, so gc has no point in placed_cell.
        """
        new_pos = []
        for index in range(len(gc.pos)):
            new_pos.append(self._point_translation(gc, index, placed_cell))
            
        return gc.__class__(gc._chord, new_pos)
    
    # sCN: places the other end of a chord correctly
    # this translates a given chord when a point within the chord was placed.
    def _gridded_chord_translation_with_point(
        self, gc: GriddedChord, point_index: int
    ) -> GriddedChord:
        """
        Return the stretched gridded chord diagram obtained when the point at
        point_index in gc is placed.
        """
        chord_placed = gc.patt[point_index] # gets the chord number of the chord being placed

        new_pos = [
            self._point_translation(gc, i, (point_index, gc.patt[point_index]))
            if gc.patt[i] != chord_placed # place the point normally if it is not in the chord being placed
            else self._point_placed_cell(gc.pos[i], i != point_index, i > point_index)
            for i in range(len(gc.pos))
        ]
        return gc.__class__(gc._chord, new_pos)

    # sCN: added boolean to tell if the placed cell is the intended one or the obligatory other endpoint
    def _point_placed_cell(self, cell: Cell, is_other_end: bool = False, is_end_sink: bool = False) -> Cell:
        """
        Return where cell gets shifted to if a point is placed in cell.

        If placed on its own row, then the y coordinate is shifted by 1.
        If placed on its own column, then the x coordinate is shifted by 1.

        If is_other_end is true, returns the cell where the point of the other 
        end of a placed chord will be added in the placed tiling
        """
        x, y = cell
        if self.own_col:
            if is_other_end:
                if is_end_sink:
                    x = x + 2
            else:
                x = x + 1

        return (x, y + 1 if self.own_row else y)
    
    # sToDo: this should probably get moved and tested in GriddedChord
    def directionmost_point(self, gc: GriddedChord, dir: int) -> Cell:
        """Returns the point (idx, val) that is the directionmost point in gc"""
        if dir == DIR_SOUTH or dir == DIR_WEST:
            return (0, 0)
        if dir == DIR_NORTH or dir == DIR_EAST:
            chord_ends = gc.chord_dict.get(len(gc) - 1)
            return (chord_ends[0], len(gc) - 1)

    # tested! -> should this be focused on all the ways a gc can be stretched if a CHORD is placed instead of a point?
    def get_multiplexes_of_chord(self, gc: GriddedChord, cell: Cell):
        """Retruns all the ways gc could be stretched if a point is placed in cell"""
        min_idx, max_idx, min_val, max_val = gc.get_bounding_box(cell)
        #print(gc.get_bounding_box(cell))

        # if not placing on own_col or own_row, we don't want to loop through multiplexes around thoses indices.
        if not self.own_col:
            max_idx = min_idx
        if not self.own_row:
            max_val = min_val
        
        multiplexes = []
        for idx in range(min_idx, max_idx + 1):
            # find the multiplexes with gc streched around placed cell
            for val in range(min_val, max_val + 1):
                multiplexes.append(self._gridded_chord_translation(gc, (idx, val)))
            
        # find point multiplexes with gc stretched with a point in placed cell
        for idx in gc.points_in_cell(cell):
            multiplexes.append(self._gridded_chord_translation_with_point(gc, idx))

        return multiplexes

    # tested!
    def get_multiplexes_of_chords(self, gcs: Iterable[GriddedChord], cell: Cell) -> List[GriddedChord]:
        """
        Return all stretched gridded chord diagrams for an iterable of gridded
        chord diagrams, assuming a point is placed in the given cell.
        """
        return list(
            chain.from_iterable(self.get_multiplexes_of_chord(gc, cell) for gc in gcs)
        )
    
    # tested!
    def added_obs(self, cell_placed: Cell, cell_end: Cell, is_end_sink: bool = False) -> Iterable[GriddedChord]:
        """Returns the obstructions that need to be added to ensure that the placed
        point is isolated in its own cell, row and column."""
        cell_x, cell_y = self._point_placed_cell(cell_placed)
        cell_end_x, _ = self._point_placed_cell(cell_end, True, is_end_sink)
        tiling_length, tiling_height = self._tiling.dimensions
        new_obs = []
        
        for x in range(tiling_length + 2):
            if (x == cell_x or x == cell_end_x):
                if not is_end_sink:
                    # obstructions to ensure point in cell is isolated -- we only want these if we are placing the second point
                    new_obs += [GriddedChord(Chord((0, 1)), ((x, cell_y), (x, cell_y))),
                                    GriddedChord(Chord((1, 0)), ((x, cell_y), (x, cell_y))),
                                    GriddedChord(Chord((0, 0)),((x, cell_y), (x, cell_y)))]
            else:
                # obstructions to ensure the place point is isolated in its own row
                new_obs.append(GriddedChord(Chord((0,)), ((x, cell_y),)))

        for y in range(tiling_height + 2):
            # we only want to add point obstructions for cells that are not the placed cell
            if y != cell_y:
                new_obs.append(GriddedChord(Chord((0,)), ((cell_x, y),)))
        #print("new obs in added obs", new_obs)
        return new_obs
    
    
    # tested!
    # sToDo: some obstructions being added currently may be redundant
    def stretched_obs(self, cell_placed: Cell) -> Iterable[GriddedChord]:
        """Returns the existing obstructions stretched over the cell that had a point placed.
        I.e., the mulitplexes of the obstructions over the cell placed"""
        return self.get_multiplexes_of_chords(self._tiling.obstructions, cell_placed)
    
    # sToDo: is_other_end might be better if changed to is_source, since we would want to place source than sink,
    # not most extreme point, other end of chord.
    # tested!
    def point_multiplex_obs(self, req_placed: GriddedChord, placed_idx: int) -> Iterable[GriddedChord]:
        """Returns the chords that need to be forbidden to ensure
        the placed point is the directionmost point for direction dir
        
        place_source says whether the source of the chord is being placed (if not it is the sink)"""
        placed_cell = req_placed.pos[placed_idx]
            
        obs = [self._gridded_chord_translation_with_point(req_placed, idx) # this is basically a point multiplex
                for idx in req_placed.points_in_cell(placed_cell)
                if idx != placed_idx] # we obviously don't want to forbid the placement of the point
        return obs
    
    # tested!
    def point_dir_obs(self, cell_placed: Cell, dir: int) -> Iterable[GriddedChord]:
        """Retruns the point obstructions to make sure there are no active cells 
        farther in the direction dir than the cell placed"""
        cell_location = self._point_placed_cell(cell_placed)
        
        obs = []

        if dir == DIR_WEST or dir == DIR_SOUTH:
            col = cell_location[0] - 1
            row = cell_location[1] - 1             

        if dir == DIR_NORTH or dir == DIR_EAST:
            row = cell_location[1] + 1

        if dir == DIR_WEST or dir == DIR_SOUTH:
            extra = 0
            if self.own_row: 
                extra += 2
            for y in range(self._tiling.dimensions[1] + extra):
                obs.append(GriddedChord(Chord((0,)), ((col, y),)))
    
        extra = 0
        if self.own_col:
            extra += 2
        for x in range(self._tiling.dimensions[0] + extra):
            obs.append(GriddedChord(Chord((0,)), ((x, row),)))        
            
        return obs

    # tested!
    def added_reqs(self, req_placed: GriddedChord, idx_placed: int, idx_end: int, is_end_sink: bool) -> Iterable[Iterable[GriddedChord]]:
        """Returns the reqs that need to be added to place at idx_placed"""
        reqs = [[GriddedChord(Chord((0,)), (self._point_placed_cell(req_placed.pos[idx_placed]),))],
                [GriddedChord(Chord((0,)), (self._point_placed_cell(req_placed.pos[idx_end], True, is_end_sink),))],
                [self._gridded_chord_translation_with_point(req_placed, idx_placed)]]
        return reqs
    
    # printed and it works
    def stretched_reqs(self, cell_placed: Cell, req_placed: GriddedChord) -> list[list[GriddedChord]]:
        stretched_req_lists = []
        reqs = list(self._tiling.requirements)
        reqs.remove((req_placed,))
        for req_list in reqs:
            stretched_req_lists.append(self.get_multiplexes_of_chords(req_list, cell_placed))

        return stretched_req_lists
    
    def stretch_links(self, cell_placed: Cell) -> list[list[Cell]]:
        new_links = []
        for link in self._tiling.linkages:
            assert link
            link_xs = [cell[0] for cell in link]
            link_ys = [cell[1] for cell in link]

            min_x = min(link_xs)
            min_y = min(link_ys)
            max_x = max(link_xs)
            max_y = max(link_ys)

            if min_x > cell_placed[0] and self.own_col:
                min_x += 2
            if max_x > cell_placed[0] and self.own_col:
                max_x += 2

            if min_y > cell_placed[0] and self.own_row:
                min_y += 2
            if max_y > cell_placed[0] and self.own_row:
                max_y += 2

            new_link = product(range(min_x, max_x + 1), range(min_y, max_y + 1))
            new_links.append(new_link)

        return new_links

    # sToDo: i wonder if there is double calculation going on. I.e., if some of the variables calculated
    # here are calculated again in methods. This should be double checked and reduced.
    def place_point(self, req_placed: GriddedChord, dir: int, place_source: bool) -> Tiling:

        _, val_placed = self.directionmost_point(req_placed, dir)

        if place_source:
            placed_idx, other_idx = req_placed.chord_dict[val_placed]
        else: 
            other_idx, placed_idx = req_placed.chord_dict[val_placed]
        
        placed_cell = req_placed.pos[placed_idx]
        end_cell = req_placed.pos[other_idx]

        # Note that if place_source is true, we are placing a source chord, so the other end of the chord is a sink
        obs = []
        obs += self.added_obs(placed_cell, end_cell, place_source)

        
        # we don't want to add the direction obstructions if we aren't placing the source
        if place_source: 
            obs += self.point_dir_obs(placed_cell, dir)

        obs += self.stretched_obs(placed_cell)

        obs += self.point_multiplex_obs(req_placed, placed_idx)

        reqs = []
        reqs += self.added_reqs(req_placed, placed_idx, other_idx, place_source)
        reqs += self.stretched_reqs(placed_cell, req_placed)

        #links = []
        #links += self.stretch_links(placed_cell)

        return Tiling(obs, reqs, tuple(), self._tiling.assumptions, remove_empty_rows_and_cols=False)
    
    def place_chord(self, req_placed: GriddedChord, dir: int) -> Tiling:
        source_placed_tiling = self.place_point(req_placed, dir, True)

        source_placement_class = self.__class__(source_placed_tiling, False, self.own_col, self.directions)
        _, val_placed = source_placement_class.directionmost_point(req_placed, dir)
        place_source_idx, place_sink_idx = req_placed.chord_dict[val_placed]

        source_placed_req = self._gridded_chord_translation_with_point(req_placed, place_source_idx)

        chord_placed = source_placement_class.place_point(source_placed_req, dir, False)
        chord_placed._simplify()
        chord_placed._remove_empty_rows_and_cols()
        return chord_placed
    
    def place_chords(self, req_list_to_place: tuple[GriddedChord], dir: int):
        res = []
        for req in req_list_to_place:
            res.append(self.place_chord(req, dir))
        return tuple(res)


non_crossing = Tiling(
        obstructions=(
            GriddedChord(Chord((0,1,0,1)), ((0, 0),) * 4),
        ),
        requirements=(
            [GriddedChord(Chord((0,0)), ((0,0), (0,0)))],
        )
    )    
'''placement_nc = RequirementPlacement(non_crossing)
print()
place = placement_nc.place_chord(GriddedChord(Chord((0,0)), ((0,0), (0,0))), 3)
place._simplify()
place._remove_empty_rows_and_cols()
already_placed = RequirementPlacement(place)
print(already_placed)
print(already_placed.chord_already_placed(GriddedChord(Chord((0,0)), ((0,0), (2,0))), 3))'''
