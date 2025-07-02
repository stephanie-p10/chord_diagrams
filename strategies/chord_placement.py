import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from typing import (Any, Callable, Dict, FrozenSet, Iterable, Iterator, List, Optional, Set, Tuple)

#import tilings.strategies as strat
from comb_spec_searcher import DisjointUnionStrategy, StrategyFactory
from comb_spec_searcher.exception import StrategyDoesNotApply
from comb_spec_searcher.strategies import Rule
from comb_spec_searcher.strategies.strategy import VerificationStrategy
from chords import Chord, GriddedChord
from tiling import Tiling

Cell = Tuple[int, int]
ListRequirement = Tuple[GriddedChord, ...]

'''
# copy and pasted code that might be useful later when more functionality is implements (like RowColMap)
class RequirementPlacementStrategy(DisjointUnionStrategy[Tiling, GriddedChord]):
    def __init__(
        self,
        gps: Iterable[GriddedChord],
        indices: Iterable[int],
        direction: int,
        own_col: bool = True,
        own_row: bool = True,
        ignore_parent: bool = False,
        include_empty: bool = False,
    ):
        self.gps = tuple(gps)
        self.indices = tuple(indices)
        self.direction = direction
        self.own_row, self.own_col = own_row, own_col
        self.include_empty = include_empty
        self._placed_cells = tuple(
            sorted(set(gp.pos[idx] for idx, gp in zip(self.indices, self.gps)))
        )
        possibly_empty = self.include_empty or len(self.gps) > 1
        super().__init__(ignore_parent=ignore_parent, possibly_empty=possibly_empty)

    def _placed_cell(self, idx: int) -> Cell:
        """Return the cell placed given the index of the child."""
        return self._placed_cells[idx]

    def _child_idx(self, idx: int):
        """Return the index of the child given the index of gps placed into."""
        return self._placed_cells.index(self.gps[idx].pos[self.indices[idx]])

    def placement_class(self, tiling: Tiling) -> RequirementPlacement:
        return RequirementPlacement(tiling, own_col=self.own_col, own_row=self.own_row)

    def decomposition_function(self, comb_class: Tiling) -> Tuple[Tiling, ...]:
        placement_class = self.placement_class(comb_class)
        placed_tilings = placement_class.place_point_of_req(
            self.gps, self.indices, self.direction
        )
        if self.include_empty:
            return (comb_class.add_obstructions(self.gps),) + placed_tilings
        return placed_tilings

    def extra_parameters(
        self, comb_class: Tiling, children: Optional[Tuple[Tiling, ...]] = None
    ) -> Tuple[Dict[str, str], ...]:
        if not comb_class.extra_parameters:
            return super().extra_parameters(comb_class, children)
        if children is None:
            children = self.decomposition_function(comb_class)
            if children is None:
                raise StrategyDoesNotApply("Strategy does not apply")
        algo = self.placement_class(comb_class)
        extra_parameters: Tuple[Dict[str, str], ...] = tuple({} for _ in children)
        if self.include_empty:
            child = children[0]
            for assumption in comb_class.assumptions:
                mapped_assumption = child.forward_map.map_assumption(
                    assumption
                ).avoiding(child.obstructions)
                if mapped_assumption.gps:
                    parent_var = comb_class.get_assumption_parameter(assumption)
                    child_var = child.get_assumption_parameter(mapped_assumption)
                    extra_parameters[0][parent_var] = child_var
        for idx, (cell, child) in enumerate(
            zip(self._placed_cells, children[1:] if self.include_empty else children)
        ):
            mapped_assumptions = [
                child.forward_map.map_assumption(ass).avoiding(child.obstructions)
                for ass in algo.stretched_assumptions(cell)
            ]
            for assumption, mapped_assumption in zip(
                comb_class.assumptions, mapped_assumptions
            ):
                if mapped_assumption.gps:
                    parent_var = comb_class.get_assumption_parameter(assumption)
                    child_var = child.get_assumption_parameter(mapped_assumption)
                    extra_parameters[idx + 1 if self.include_empty else idx][
                        parent_var
                    ] = child_var
        return extra_parameters

    def direction_string(self):
        if self.direction == DIR_EAST:
            return "rightmost"
        if self.direction == DIR_NORTH:
            return "topmost"
        if self.direction == DIR_WEST:
            return "leftmost"
        if self.direction == DIR_SOUTH:
            return "bottommost"

    def formal_step(self):
        placing = f"placing the {self.direction_string()} "
        if not (self.own_row and self.own_col):
            placing = f"partially {placing}"
        if len(self.gps) == 1:
            gp = self.gps[0]
            index = self.indices[0]
            if len(gp) == 1:
                return placing + f"point in cell {gp.pos[index]}"
            if gp.is_localized():
                return (
                    f"{placing}{(index, gp.patt[index])} point in "
                    f"{gp.patt} in cell {gp.pos[index]}"
                )
            return f"{placing}{(index, gp.patt[index])} point in {gp}"
        if all(len(gp) == 1 for gp in self.gps):
            col_indices = set(x for x, _ in [gp.pos[0] for gp in self.gps])
            if len(col_indices) == 1:
                return f"{placing}point in column {col_indices.pop()}"
            row_indices = set(y for _, y in [gp.pos[0] for gp in self.gps])
            if len(row_indices) == 1:
                return f"{placing}point in row {row_indices.pop()}"
        return (
            f"{placing}point at indices {self.indices} from the requirement "
            f"({', '.join(map(str, self.gps))})"
        )

    def backward_cell_map(self, placed_cell: Cell, cell: Cell) -> Cell:
        x, y = cell
        if self.own_col and x > placed_cell[0] + 1:
            x -= 2
        elif self.own_col and x == placed_cell[0] + 1:
            x -= 1
        if self.own_row and y > placed_cell[1] + 1:
            y -= 2
        elif self.own_row and y == placed_cell[1] + 1:
            y -= 1
        return x, y

    def forward_gp_map(self, gp: GriddedChord, forced_index: int) -> GriddedChord:
        new_pos: List[Cell] = []
        forced_val = gp.patt[forced_index]
        for idx, (x, y) in enumerate(gp.pos):
            if idx == forced_index:
                if self.own_col:
                    x += 1
                if self.own_row:
                    y += 1
                new_pos.append((x, y))
            else:
                val = gp.patt[idx]
                if self.own_col and idx >= forced_index:
                    x += 2
                if self.own_row and val >= forced_val:
                    y += 2
                new_pos.append((x, y))
        return GriddedChord(gp.patt, new_pos)

    def backward_map(
        self,
        comb_class: Tiling,
        objs: Tuple[Optional[GriddedChord], ...],
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Iterator[GriddedChord]:
        if children is None:
            children = self.decomposition_function(comb_class)
        idx = DisjointUnionStrategy.backward_map_index(objs)
        gp: GriddedChord = children[idx].backward_map.map_gp(
            cast(GriddedChord, objs[idx])
        )
        if self.include_empty:
            if idx == 0:
                yield gp
                return
            idx -= 1
        placed_cell = self._placed_cell(idx)
        yield GriddedChord(
            gp.patt, [self.backward_cell_map(placed_cell, cell) for cell in gp.pos]
        )

    def forward_map(
        self,
        comb_class: Tiling,
        obj: GriddedChord,
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Tuple[Optional[GriddedChord], ...]:
        indices = obj.forced_point_of_requirement(
            self.gps, self.indices, self.direction
        )
        if children is None:
            children = self.decomposition_function(comb_class)
        if indices is None:
            return (children[0].forward_map.map_gp(obj),) + tuple(
                None for _ in range(len(children) - 1)
            )
        gps_index, forced_index = indices
        child_index = self._child_idx(gps_index)
        if self.include_empty:
            child_index += 1
        gp = self.forward_gp_map(obj, forced_index)
        return (
            tuple(None for _ in range(child_index))
            + (children[child_index].forward_map.map_gp(gp),)
            + tuple(None for _ in range(len(children) - 1))
        )

    def __str__(self) -> str:
        return "requirement placement strategy"

    def __repr__(self) -> str:
        return (
            f"RequirementPlacementStrategy(gps={self.gps}, "
            f"indices={self.indices}, direction={self.direction}, "
            f"own_col={self.own_col}, own_row={self.own_row}, "
            f"ignore_parent={self.ignore_parent}, "
            f"include_empty={self.include_empty})"
        )

    def to_jsonable(self) -> dict:
        """Return a dictionary form of the strategy."""
        d: dict = super().to_jsonable()
        d.pop("workable")
        d.pop("inferrable")
        d.pop("possibly_empty")
        d["gps"] = tuple(gp.to_jsonable() for gp in self.gps)
        d["indices"] = self.indices
        d["direction"] = self.direction
        d["own_col"] = self.own_col
        d["own_row"] = self.own_row
        d["include_empty"] = self.include_empty
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "RequirementPlacementStrategy":
        gps = tuple(GriddedChord.from_dict(gp) for gp in d.pop("gps"))
        return cls(gps=gps, **d)'''

class RequirementPlacementStrategy(DisjointUnionStrategy[Tiling, GriddedChord]):
    def __init__(
        self,
        gps: Iterable[GriddedChord],
        indices: Iterable[int],
        direction: int,
        own_col: bool = True,
        own_row: bool = True,
        ignore_parent: bool = False,
        include_empty: bool = False,
    ):
        self.gps = tuple(gps)
        self.indices = tuple(indices)
        self.direction = direction
        self.own_row, self.own_col = own_row, own_col
        self.include_empty = include_empty
        self._placed_cells = tuple(
            sorted(set(gp.pos[idx] for idx, gp in zip(self.indices, self.gps)))
        )
        possibly_empty = self.include_empty or len(self.gps) > 1
        super().__init__(ignore_parent=ignore_parent, possibly_empty=possibly_empty)