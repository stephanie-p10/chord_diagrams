"""Algorithms for inserting/placing requirements into a tiling.

This module implements the mechanics behind requirement placement strategies:
given a tiling and a target gridded chord (often representing a “single point”
or “single chord” requirement), compute the derived obstructions/requirements
and linkage updates induced by choosing a particular placement.

The code is cache-heavy because strategies may query many candidate placements
repeatedly; the caches are keyed by cells and (requirement, index) pairs.
"""

import sys
from pathlib import Path
 
_src_root = Path(__file__).resolve().parents[2]  # .../src
if str(_src_root) not in sys.path:
    sys.path.insert(0, str(_src_root))
 
from src.chord_diagrams.misc import DIR_EAST, DIR_NONE, DIR_NORTH, DIR_SOUTH, DIR_WEST, DIRS
from src.chord_diagrams.assumptions import TrackingAssumption
from src.chord_diagrams.chords import Chord, GriddedChord
from src.chord_diagrams.tiling import Tiling


from itertools import chain, filterfalse, product, chain
from typing import TYPE_CHECKING, Dict, FrozenSet, Iterable, List, Optional, Tuple
from collections import Counter

#from .simplify import SimplifyObstructionsAndRequirements

Cell = Tuple[int, int]
Dir = int
ListRequirement = List[GriddedChord]
Linkages = List[Cell]
LinkCache = Dict[Cell, List[Linkages]]
ObsCache = Dict[Cell, List[GriddedChord]]
ReqsCache = Dict[Cell, List[ListRequirement]]
AssumpCache = Dict[Cell, List[TrackingAssumption]]

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
    
    def _point_translation(self, gc: GriddedChord, idx: int, point_placed: Cell) -> Cell:
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
    
    def new_empty_cells(self, cell_placed: Cell, cell_end: Cell, dir: int, is_end_sink: bool = False) -> Iterable[Cell]:
        cell_x, cell_y = self._point_placed_cell(cell_placed)
        cell_end_x, _ = self._point_placed_cell(cell_end, True, is_end_sink)
        tiling_length, tiling_height = self._tiling.dimensions
        empty_cells = []

        # cells in same row of placed chord that don't contain the placed chord must be empty
        for x in range(tiling_length + 2):
            if not (x == cell_x or x == cell_end_x):
                empty_cells.append((x, cell_y),)

        # cells that are in the same col of the placed point must be empty
        for y in range(tiling_height + 2):
            if y != cell_y:
                empty_cells.append((cell_x, y),)

        # if the other end is a sink, then this point is a source, and we want to place direction obs
        if is_end_sink:
            """Retruns the point obstructions to make sure there are no active cells 
            farther in the direction dir than the cell placed"""
            cell_location = self._point_placed_cell(cell_placed)

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
                    empty_cells.append((col, y),)
        
            extra = 0
            if self.own_col:
                extra += 2
            for x in range(self._tiling.dimensions[0] + extra):
                empty_cells.append((x, row),)      

        return empty_cells
    
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
                new_obs.append(GriddedChord(Chord((0,)), ((x, cell_y),))) # added to empty cells

        for y in range(tiling_height + 2):
            # we only want to add point obstructions for cells that are not the placed cell
            if y != cell_y:
                new_obs.append(GriddedChord(Chord((0,)), ((cell_x, y),))) # added to empty cells
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
                obs.append(GriddedChord(Chord((0,)), ((col, y),))) # added to empty cells
    
        extra = 0
        if self.own_col:
            extra += 2
        for x in range(self._tiling.dimensions[0] + extra):
            obs.append(GriddedChord(Chord((0,)), ((x, row),)))  # added to empty cells   
            
        return obs

    # tested!
    def added_reqs(self, req_placed: GriddedChord, idx_placed: int, idx_end: int, is_end_sink: bool) -> Iterable[Iterable[GriddedChord]]:
        """Returns the reqs that need to be added to place at idx_placed"""
        reqs = [[GriddedChord(Chord((0,)), (self._point_placed_cell(req_placed.pos[idx_placed]),))],
                [GriddedChord(Chord((0,)), (self._point_placed_cell(req_placed.pos[idx_end], True, is_end_sink),))],
                [self._gridded_chord_translation_with_point(req_placed, idx_placed)]]
        return reqs
    
    # printed and it works
    def stretched_reqs(self, cell_placed: Cell, req_placed: GriddedChord) -> List[List[GriddedChord]]:
        stretched_req_lists = []
        reqs = list(self._tiling.requirements)
        reqs.remove((req_placed,))
        for req_list in reqs:
            stretched_req_lists.append(self.get_multiplexes_of_chords(req_list, cell_placed))

        return stretched_req_lists
    
    def stretch_links(self, cell_placed: Cell) -> List[List[Cell]]:
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

        return Tiling(obs, reqs, tuple(), self._tiling.assumptions, remove_empty_rows_and_cols=False) #why is remove_empty_rc false?
    
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
    
    def place_chords(self, req_list_to_place: Tuple[GriddedChord, ...], dir: int):
        res = []
        for req in req_list_to_place:
            res.append(self.place_chord(req, dir))
        return tuple(res)


class ChordPlacement:
    """
    Container class for placing a single chord in a singleton requirment

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

    def __init__(
        self,
        tiling: "Tiling",
        own_row: bool = True,
        own_col: bool = True,
        dirs: Iterable[int] = tuple(DIRS),
    ):
        self._tiling = tiling
        self.own_row = own_row
        self.own_col = own_col
        
        # Checks that parameters given give a valid placement
        if not own_row and not own_col:
            raise ValueError("Must place on own row or on own column.")

        assert all(d in DIRS for d in dirs), "Got an invalid direction"

        if self.own_row and self.own_col:
            self.directions = frozenset(DIRS)
        elif self.own_row:
            self.directions = frozenset((DIR_NORTH, DIR_SOUTH))
        elif self.own_col:
            self.directions = frozenset((DIR_EAST, DIR_WEST))
        self.directions = frozenset(dirs).intersection(self.directions)
        assert self.directions, "No direction to place"

    
    
    def chord_already_placed(self, gc: GriddedChord, chord: int) -> bool:
        """
        Determines if chord is already isolated in its own row and own column on gc
        """
        source_cell, sink_cell = gc.chord_dict[chord]

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

    def chords_already_placed(self, gcs: Iterable[GriddedChord], chords: Iterable[int]) -> bool:
        all_placed = all(self.chord_already_placed(gc, chord) for gc, chord in zip(gcs, chords))
        return all_placed
    
    def _point_translation(
            self, gc: GriddedChord, idx: int, placed_cell: Cell
            ) -> Tuple[int, int]:
        """Finds the cell that the point in gc at idx gets shifted to if a 
        point was placed in placed_cell.
        
        Note that placed_cell is described in terms of a position in gc; 
        Eg. if placed_cell = (x, y), a cell in gc gets shifted around a point
        being placed between indices x-1 and x in gc, and chord values y-1 and y"""

        x, y = gc.pos[idx]
        x = x + 2 if self.own_col and idx >= placed_cell[0] else x
        y = y + 2 if (self.own_row and gc.patt[idx] >= placed_cell[1]) else y
        return (x, y)

    def _gridded_chord_translation(
        self, gc: GriddedChord, placed_cell: Cell
    ) -> GriddedChord:
        """
        Return the gridded chord diagram with its positions translated around
        a point placed at placed_cell = (idx, val). Assumes no cell in gc was 
        the cell being placed, so gc has no point in placed_cell.

        Note that placed_cell descibes where a point is being placed in terms
        of gc; Eg. if placed_cell = (x,y), a point is being placed between index 
        x-1 and x and between chord value y-1 and y. 
        """
        new_pos = []
        for index in range(len(gc.pos)):
            new_pos.append(self._point_translation(gc, index, placed_cell))
            
        return gc.__class__(gc._chord, new_pos)
    
    def _chord_placed_cell(self, cell: Cell, is_other_end: bool = False, is_end_sink: bool = False) -> Cell:
        """
        Return where cell gets shifted to when a chord is placed at cell[1].

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
    
    def _gridded_chord_translation_with_point(
        self, gc: GriddedChord, point_index: int
    ) -> GriddedChord:
        """
        Return the resulting gridded chord diagram from isolating the point
        at point_index.

        These are point multiplexes around the point point_index in gc. 
        """
        chord_placed = gc.patt[point_index] # gets the chord number of the chord being placed

        new_pos = [
            self._point_translation(gc, i, (point_index, gc.patt[point_index]))
            if gc.patt[i] != chord_placed # place the point normally if it is not in the chord being placed
            else self._chord_placed_cell(gc.pos[i], i != point_index, i > point_index)
            for i in range(len(gc.pos))
        ]

        return gc.__class__(gc._chord, new_pos)
    
    def get_multiplexes_of_gc(self, gc: GriddedChord, cell: Cell):
        """Returns all the ways gc could be stretched if a point is placed in cell

        This is the set M_{cell}(gc) as in ABCNPU. However, note that only multiplexes
        with at most one point in the placed cell are added, since the others would
        be obstructed when cell is required to be a point cell."""

        min_idx, max_idx, min_val, max_val = gc.get_bounding_box(cell)

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
    
    def get_multiplexes_of_gcs(self, gcs: Iterable[GriddedChord], cell: Cell) -> List[GriddedChord]:
        """
        Return all stretched gridded chord diagrams for an iterable of gridded
        chord diagrams, assuming a point is placed in the given cell.

        This is the set M_{cell}(gcs) in the notation of ABCNPU. 
        """
        return list(
            chain.from_iterable(self.get_multiplexes_of_gc(gc, cell) for gc in gcs)
        )
    
    def is_directionmost(self, dir: int, a: Cell, b: Cell) -> bool:
        """
        Returns if cell a is further in direction dir then cell b"""
        if dir == DIR_NORTH:
            return a[1] > b[1]
        if dir == DIR_EAST:
            return a[0] > b[0]
        if dir == DIR_SOUTH:
            return a[1] < b[1]
        if dir == DIR_WEST:
            return a[0] < b[0]
        raise ValueError("must give valid direction")
    
    # NOT TESTED TODO
    def translate_linkage(self, cell_placed: Cell, linkage: List[Cell]):
        new_linkage = []
        for cell in linkage:
            expand_x = 0
            expand_y = 0
            if cell_placed[0] == cell[0]:
                expand_x = 2
            if cell_placed[1] == cell[1]:
                expand_y = 2

            corrected_x = cell[0]
            if cell_placed[0] < cell[0]:
                corrected_x += 2

            corrected_y = cell[1]
            if cell_placed[1] < cell[1]:
                corrected_y += 2

            for i in range(corrected_x, corrected_x + expand_x + 1):
                for j in range(corrected_y, corrected_y + expand_y + 1):
                    new_linkage.append((i, j))

    def _isolate_point_obs(self, isolated_cell: Cell, end_cell: Cell, dimensions: Tuple[int, int]):
        isolate_chord = [GriddedChord(Chord((0, 0)), (isolated_cell, isolated_cell)),
                       GriddedChord(Chord((1, 0)), (isolated_cell, isolated_cell)),
                       GriddedChord(Chord((0, 1)), (isolated_cell, isolated_cell)),
                       
                       GriddedChord(Chord((0, 0)), (end_cell, end_cell)),
                       GriddedChord(Chord((1, 0)), (end_cell, end_cell)),
                       GriddedChord(Chord((0, 1)), (end_cell, end_cell)),] 
        
        width, height = dimensions

        isolate_col = [GriddedChord(Chord((0,)), ((isolated_cell[0], y),)) for y in range(height + 2)]
        isolate_row = [GriddedChord(Chord((0,)), ((x, isolated_cell[1]),)) for x in range(width + 2)]
        
        isolate_col.remove(GriddedChord(Chord((0,)), (isolated_cell,)))
        isolate_row.remove(GriddedChord(Chord((0,)), (isolated_cell,)))
        isolate_row.remove(GriddedChord(Chord((0,)), (end_cell,)))

        return list(chain(isolate_chord, isolate_row, isolate_col))


    def place_point(self, req_placed: GriddedChord, idx_to_place: int, idx_other_end: int, dir: int, tiling_placed: Tiling = None, remove_empty: bool = True):
        """Places the point at idx in req_placed in the direction dir on tiling_placed."""
        if tiling_placed == None:
            tiling_placed = self._tiling
        #print("tiling to place", tiling_placed)
        cell_placed = req_placed.pos[idx_to_place]  # 'c' in notes

        translated_req_placed = self._gridded_chord_translation_with_point(req_placed, idx_to_place)  # {m_c^I(r)} in notes
        translated_placed_cell = translated_req_placed.pos[idx_to_place]
        translated_end_cell = translated_req_placed.pos[idx_other_end]

        isolate_chord_obs = self._isolate_point_obs(translated_placed_cell, translated_end_cell, tiling_placed.dimensions) # A \cup B_1 \cup B_2 in notes

        obs_multiplexes = self.get_multiplexes_of_gcs(tiling_placed.obstructions, cell_placed)  # O' := M_c(O) in notes

        multiplexes: List[GriddedChord] = self.get_multiplexes_of_gc(req_placed, cell_placed) # S in notes
        direction_obs = []
        for multiplex in multiplexes:
            if self.is_directionmost(dir, multiplex.pos[idx_to_place], translated_placed_cell):
                direction_obs.append(multiplex)

        translated_req_lists = [] # {R_2',...,R_k'} in notes.
        # However, the requirement list being placed may not be the first one WLOG like on paper
        for req_list in tiling_placed.requirements:
            if not (len(req_list) == 1 and req_placed in req_list):
                translated_req_lists.append(self.get_multiplexes_of_gcs(req_list, cell_placed))

        
        # makes sure there is actually a requirement that was placed. 
        assert len(translated_req_lists) < len(tiling_placed.requirements)

        single_point_req = GriddedChord(Chord((0,)), (translated_req_placed.pos[idx_to_place],)) # C in notes

        new_linkages = [self.translate_linkage(cell_placed, linkage) for linkage in tiling_placed.linkages]

        new_obs = list(chain(obs_multiplexes, isolate_chord_obs, direction_obs))

        new_reqs = list(chain([(translated_req_placed,)], translated_req_lists, [(single_point_req,)]))

        return Tiling(new_obs, new_reqs, new_linkages, remove_empty_rows_and_cols = remove_empty)
    
    def place_chord(self, req_placed: GriddedChord, chord_to_place: int, dir: int):
        source_idx, sink_idx = req_placed.chord_dict[chord_to_place]
        idx_place_first = source_idx
        idx_place_second = sink_idx
        if dir == DIR_EAST or dir == DIR_WEST:
            idx_place_first = sink_idx
            idx_place_second = source_idx

        # Tiling after the first point of chord_to_place has been placed. Note that empty rows and columns remain 
        # to make it easier to calculate the updated position of the requirement being placed. 
        placed_first_point = self.place_point(req_placed, idx_place_first, idx_place_second, dir, remove_empty=False)

        # This is the requirement being placed updated to the position it would be after the first point has been placed.
        translated_req_placed = self._gridded_chord_translation_with_point(req_placed, idx_place_first) 

        # Tiling after the second point of the chord being placed has been placed.
        placed_second_point = self.place_point(translated_req_placed, idx_place_second, idx_place_first, dir, placed_first_point, remove_empty=True)

        return placed_second_point

    def place_chords(self, req_list_to_place: Tuple[GriddedChord, ...], chords_to_place: Iterable[int], dir: int):
        res = []
        for req, chord in zip(req_list_to_place, chords_to_place):
            res.append(self.place_chord(req, chord, dir))
        return tuple(res)


# first pattern in theorem 3.1.2 class
c1 = Chord((0, 1, 2, 0, 1, 2))
# second pattern in theorem 3.1.2 class
c2 = Chord((0, 1, 2, 0, 2, 1))
t5 = Tiling(obstructions=(GriddedChord.single_chord(((0, 0), (0, 0))),
                                      GriddedChord.single_chord(((2, 0), (2, 0))), 
                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                      GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0))),

                                      GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                      GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                      GriddedChord(c1, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                      GriddedChord(c1, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),

                                      GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                      GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                      GriddedChord(c2, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                      GriddedChord(c2, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                
                                      GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (3, 1), (3, 1))),
                                      GriddedChord(Chord((0, 1, 1, 0)), ((1, 1), (1, 1), (3, 1), (3, 1))),),
                        requirements=((GriddedChord.single_chord(((0, 0), (2, 0))),),
                                      (GriddedChord.single_chord(((1, 1), (3, 1))),)))

t3 = Tiling(obstructions=(GriddedChord(c1, ((0, 0),) *6), GriddedChord(c2, ((0, 0),) *6)),
            requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),),))

non_crossing = Tiling(
    obstructions=(
        GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),
    ),
    requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),),),
)

#if __name__ == "__main__":
placement_nc = ChordPlacement(non_crossing)
placement = ChordPlacement(t3)

GriddedChord.single_chord([(0, 0), (0, 0)])
GriddedChord(Chord((0, 0)), ((0, 0), (0, 0)))

#print("placement:")

t = placement_nc.place_chord(GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))), 0, DIR_SOUTH)
#print(t)
t._simplify()
#print(t)
#print([t])
