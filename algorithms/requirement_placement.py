import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from itertools import chain, filterfalse
from typing import TYPE_CHECKING, Dict, FrozenSet, Iterable, List, Optional, Tuple

from permuta.misc import DIR_EAST, DIR_NONE, DIR_NORTH, DIR_SOUTH, DIR_WEST, DIRS
from assumptions import TrackingAssumption
from chords import Chord, GriddedChord
from steph_chords.tiling import Tiling

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
            - `permuta.misc.DIR_NORTH`
            - `permuta.misc.DIR_SOUTH`
            - `permuta.misc.DIR_EAST`
            - `permuta.misc.DIR_WEST`
    """

    def __init__(
        self,
        tiling: "Tiling",
        own_row: bool = True,
        own_col: bool = True,
        dirs: Iterable[int] = tuple(DIRS),
    ):
        if not own_row and not own_col:
            raise ValueError("Must place on own row or on own column.")
        assert all(d in DIRS for d in dirs), "Got an invalid direction"
        self._tiling = tiling
        self._point_row_cells = self._tiling_point_row_cells()
        self._point_col_cells = self._tiling_point_col_cells()
        self.own_row = own_row
        self.own_col = own_col
        self._stretched_obstructions_cache: ObsCache = {}
        self._stretched_requirements_cache: ReqsCache = {}
        self._stretched_assumptions_cache: AssumpCache = {}
        self._stretched_linkages_cache: Linkages
        if self.own_row and self.own_col:
            self.directions = frozenset(DIRS)
        elif self.own_row:
            self.directions = frozenset((DIR_NORTH, DIR_SOUTH))
        elif self.own_col:
            self.directions = frozenset((DIR_EAST, DIR_WEST))
        self.directions = frozenset(dirs).intersection(self.directions)
        assert self.directions, "No direction to place"

    def _tiling_point_col_cells(self) -> FrozenSet[Cell]:
        """
        The point cells of the tilings that are the only active cell in their
        column.
        """
        return frozenset(
            filter(self._tiling.only_cell_in_col, self._tiling.point_cells) # filters point_cells by if they are the only cells in their col
        )

    def _tiling_point_row_cells(self) -> FrozenSet[Cell]:
        """
        The point cells of the tilings that are the only active cell in their
        row.
        """
        return frozenset(
            filter(self._tiling.only_cell_in_row, self._tiling.point_cells) # filters point_cells by if they are the only cells in their row.
        )

    '''def already_placed(
        self, gcs: Iterable[GriddedChord], indices: Iterable[int]
    ) -> bool:
        """
        Determine if this set gridded chord diagrams is already placed.
        """
        cells = set(gc.pos[idx] for gc, idx in zip(gcs, indices)) # not sure what this is doing
        if len(cells) != 1:
            return False
        cell = cells.pop()
        full_placement = self.own_row and self.own_col
        if full_placement:
            return cell in self._point_col_cells and cell in self._point_row_cells
        if self.own_row:  # Only placing in own row
            return cell in self._point_row_cells
        if self.own_col:  # Only placing in own column
            return cell in self._point_col_cells
        raise Exception("Not placing at all!!")'''
    
    def _point_translation(self, gc: GriddedChord, index: int, placed_point: Cell) -> Cell:
        """
        Return the translated position of the cell at the given index.

        The translation assumes that there has been a point placed in the
        position (i, j) = placed_point where this corresponds to the index and
        value within the pattern of the gridded chord diagram gc. 
        (Cartesian point of chord diagram, not cell in gridded chord)

        If the newly placed point is assumed to be put on the new column we
        have that the cell is expanded like:
            -      - - -
           | | -> | |o| |
            -      - - -
        meaning that indices to the right of i are shifted by 2.
        Similarly, for new rows we have
                   -
                  | |
            -      -
           | | -> |o|
            -      -
                  | |
                   -
        meaning that values above j are shifted by 2.
        """
        x, y = gc.pos[index]
        return (
            x + 2 if self.own_col and index >= placed_point[0] else x,
            y + 2 if (self.own_row and gc.patt[index] >= placed_point[1]) else y,
        )

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
    
    # this translates a given chord around a placed point that was not in this chord
    def _gridded_chord_translation(
        self, gc: GriddedChord, placed_cell: Cell
    ) -> GriddedChord:
        """
        Return the gridded chord diagram with all of the cells translated
        assuming that a point was placed at placed_cell = (idx, val)
        """
        newpos = [
            self._point_translation(gc, index, placed_cell) for index in range(len(gc.pos))
        ]
        return gc.__class__(gc.patt, newpos)

    # sCN: places the other end of a chord correctly
    # this translates a given chord when a point within the chord was placed.
    def _gridded_chord_translation_with_point(
        self, gc: GriddedChord, point_index: int
    ) -> GriddedChord:
        """
        Return the stretched gridded permutation obtained when the point at
        point_index in gc is placed.
        """
        # the following comment was inherited from permutations:
        # TODO: to prepare for intervals consider all ways of drawing a
        #       rectangle around point in cell.
        chord_placed = gc.patt[point_index]

        new_pos = [
            self._point_translation(gc, i, (point_index, gc.patt[point_index]))
            if gc.patt[i] != chord_placed # place the point normally if it is not in the chord being placed
            else self._placed_cell(gc.pos[point_index], i != point_index)
            for i in range(len(gc.pos))
        ]
        return gc.__class__(gc.patt, new_pos)

    # sCN: added boolean to tell if the placed cell is the intended one or the obligatory other endpoint
    def _placed_cell(self, cell: Cell, is_other_end: bool = False) -> Cell:
        """
        Return the cell in which the point will be added in the placed tiling.

        If placed on its own row, then the y coordinate is shifted by 1.
        If placed on its own column, then the x coordinate is shifted by 1.

        If is_other_end is true, returns the cell where the point of the other 
        end of a placed chord will be added in the placed tiling
        """
        x, y = cell
        return (x + 1 if self.own_col and not is_other_end else x, y + 1 if self.own_row else y)

    # sAsk, sToDo: what happens when both ends of chord in same cell?
    # sCN: changed to be the obstructions that make a chord a single point
    def _point_obstructions(self, cell: Cell) -> List[GriddedChord]:
        """
        Return the localised 11 and 12 obstruction required to ensure the
        newly placed point is a point.
        """
        placed_cell = self._placed_cell(cell)
        return [
            GriddedChord((0, 0), (placed_cell, placed_cell)),
            GriddedChord((0, 1), (placed_cell, placed_cell)),
        ]

    def _point_requirements(self, cell: Cell) -> List[ListRequirement]:
        """
        Return the requirement required to ensure that the newly placed point
        is a point.
        """
        placed_cell = self._placed_cell(cell)
        return [[GriddedChord((0,), (placed_cell,))]]

    def _stretch_gridded_chord(self, gc: GriddedChord, cell: Cell) -> Iterable[GriddedChord]:
        """
        Return all of the possible ways that a gridded chord diagram can be
        stretched assuming that a point is placed into the given cell.
        """
        # left off here. Trying to update the req placement algorithm to work with chords
        mindex, maxdex, minval, maxval = gc.get_bounding_box(cell)
        if not self.own_col:
            maxdex = mindex
        elif not self.own_row:
            maxval = minval
        # all the ways a gridded chord could be placed around a placed point (from) a different chord
        res = [
            self._gridded_chord_translation(gc, (idx, val))
            for idx in range(mindex, maxdex + 1)
            for val in range(minval, maxval + 1)
        ]
        # all the ways a gc can be placed with a placed point in cell
        for i in gc.points_in_cell(cell):
            res.append(self._gridded_chord_translation_with_point(gc, i))
        return res

    def _stretch_gridded_chords(
        self, gcs: Iterable[GriddedChord], cell: Cell
    ) -> List[GriddedChord]:
        """
        Return all stretched gridded permuations for an iterable of gridded
        permutations, assuming a point is placed in the given cell.
        """
        return list(
            chain.from_iterable(self._stretch_gridded_chord(gc, cell) for gc in gcs)
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
        a gridded permutation in gps, such that the point at idx is the placed
        cell.
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
        Return the tilings, where the placed point corresponds to the directionmost
        (the furtest in the given direction, ex: leftmost point) of an occurrence
        of any point idx, gc(idx) for gridded chord diagrams in gpc, and idx in indices
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
                if not any(o2 in o1 for o2 in filterfalse(o1.__eq__, forced_obs))
            ]
            reduced_obs = [o1 for o1 in obs if not any(o2 in o1 for o2 in forced_obs)]
            reduced_obs.extend(filterfalse(reduced_obs.__contains__, forced_obs))
            res.append(
                self._tiling.__class__(
                    reduced_obs,
                    reqs + [remaining_req],
                    links,
                    assumptions=assump,
                    already_minimized_obs=True,
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
        req = [GriddedChord((0,), (cell,)) for cell in self._tiling.cells_in_col(index)]
        return self.place_point_of_req(req, tuple(0 for _ in req), direction)

    def row_placement(self, index: int, direction: Dir) -> Tuple["Tiling", ...]:
        """
        Return the list corresponding the index row being placed in the given
        direction.
        """
        assert direction in (DIR_NORTH, DIR_SOUTH)
        req = [GriddedChord((0,), (cell,)) for cell in self._tiling.cells_in_row(index)]
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