import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from typing import (Any, Callable, Dict, FrozenSet, Iterable, Iterator, List, Optional, Set, Tuple, cast)

#import tilings.strategies as strat
from comb_spec_searcher import DisjointUnionStrategy, StrategyFactory
from comb_spec_searcher.exception import StrategyDoesNotApply
from comb_spec_searcher.strategies import Rule
from comb_spec_searcher.strategies.strategy import VerificationStrategy
from chords import Chord, GriddedChord
from tiling import Tiling
from algorithms.requirement_placement import RequirementPlacement
from misc import DIR_EAST, DIR_NONE, DIR_NORTH, DIR_SOUTH, DIR_WEST, DIRS

Cell = Tuple[int, int]
ListRequirement = Tuple[GriddedChord, ...]

class RequirementPlacementStrategy(DisjointUnionStrategy[Tiling, GriddedChord]):
    def __init__(
        self,
        gcs: Iterable[GriddedChord],
        direction: int,
        indices: Iterable[int] = None,
        own_col: bool = True,
        own_row: bool = True,
        ignore_parent: bool = False,
        include_empty: bool = False,
    ):
        self.gcs = tuple(gcs) # I think this is going to be a single chord for chord pl strategy
        self.direction = direction

        calculated_indices = []
        if indices == None:
            for gc in self.gcs:
                if self.direction == DIR_NORTH or self.direction == DIR_EAST:
                    idx = gc.chord_dict[len(gc) - 1][0]
                else:
                    idx = 0

                calculated_indices.append(idx)
        else:
            calculated_indices = indices

        self.indices = tuple(calculated_indices)
        self.own_row, self.own_col = own_row, own_col
        self.include_empty = include_empty
        self._placed_cells = tuple(
            sorted(set(gc.pos[idx] for idx, gc in zip(self.indices, self.gcs)))
        )
        possibly_empty = self.include_empty or len(self.gcs) > 1
        super().__init__(ignore_parent=ignore_parent, possibly_empty=possibly_empty)

    def _placed_cell(self, idx: int) -> Cell:
        """Return the cell placed given the index of the child."""
        return self._placed_cells[idx]

    def _child_idx(self, idx: int):
        """Return the index of the child given the index of gcs placed into."""
        return self._placed_cells.index(self.gcs[idx].pos[self.indices[idx]])
        # in the list of placed cells, what is the index of the cell placed in the gc at idx in gcs

    def placement_class(self, tiling: Tiling) -> RequirementPlacement:
        return RequirementPlacement(tiling, own_col=self.own_col, own_row=self.own_row)

    def decomposition_function(self, comb_class: Tiling) -> Tuple[Tiling, ...]:
        placement_class = self.placement_class(comb_class)
        # if the gcs are a valid requirement that can be placed and they have not already been placed
        if self.gcs in comb_class.requirements and not placement_class.chords_already_placed(self.gcs, self.direction):
            placed_tilings = placement_class.place_chords(self.gcs, self.direction)
            if self.include_empty:
                return (comb_class.add_obstructions(self.gcs),) + placed_tilings
            return placed_tilings
        else:
            raise StrategyDoesNotApply

    '''def extra_parameters(
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
                if mapped_assumption.gcs:
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
        return extra_parameters'''

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
        if len(self.gcs) == 1:
            gc = self.gcs[0]
            index = self.indices[0]
            if len(gc) == 1:
                return placing + f"point in cell {gc.pos[index]}"
            if gc.is_localized():
                return (
                    f"{placing}{(index, gc.patt[index])} point in "
                    f"{gc.patt} in cell {gc.pos[index]}"
                )
            return f"{placing}{(index, gc.patt[index])} point in {gc}"
        if all(len(gc) == 1 for gc in self.gcs):
            col_indices = set(x for x, _ in [gc.pos[0] for gc in self.gcs])
            if len(col_indices) == 1:
                return f"{placing}point in column {col_indices.pop()}"
            row_indices = set(y for _, y in [gc.pos[0] for gc in self.gcs])
            if len(row_indices) == 1:
                return f"{placing}point in row {row_indices.pop()}"
        return (
            f"{placing}point at indices {self.indices} from the requirement "
            f"({', '.join(map(str, self.gcs))})"
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

    def forward_gc_map(self, gc: GriddedChord, forced_index: int) -> GriddedChord:
        """Returns what gc would be mapped to in the child class where forced_index is placed"""
        new_pos: List[Cell] = []
        forced_val = gc.patt[forced_index]
        for idx, (x, y) in enumerate(gc.pos):
            if gc.patt[idx] == forced_val:
                if self.own_col and idx == forced_index:
                    x += 1
                if self.own_row:
                    y += 1
                new_pos.append((x, y))
            else:
                val = gc.patt[idx]
                if self.own_col and idx >= forced_index:
                    x += 2
                if self.own_row and val >= forced_val:
                    y += 2
                new_pos.append((x, y))
        return GriddedChord(gc.patt, new_pos)

    def backward_map(
        self,
        comb_class: Tiling,
        objs: Tuple[Optional[GriddedChord], ...],
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Iterator[GriddedChord]:
        if children is None:
            children = self.decomposition_function(comb_class)
        idx = DisjointUnionStrategy.backward_map_index(objs)
        gc: GriddedChord = children[idx].backward_map.map_gc(cast(GriddedChord, objs[idx]))
        if self.include_empty:
            if idx == 0:
                yield gc
                return
            idx -= 1
        placed_cell = self._placed_cell(idx)
        yield GriddedChord(
            gc.patt, [self.backward_cell_map(placed_cell, cell) for cell in gc.pos]
        )

    def forward_map(
        self,
        comb_class: Tiling,
        obj: GriddedChord,
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Tuple[Optional[GriddedChord], ...]:
        indices = obj.forced_point_of_requirement(
            self.gcs, self.indices, self.direction
        )
        if children is None:
            children = self.decomposition_function(comb_class)
        if indices is None:
            return (children[0].forward_map.map_gc(obj),) + tuple(
                None for _ in range(len(children) - 1)
                )
        gcs_index, forced_index = indices
        child_index = self._child_idx(gcs_index)
        if self.include_empty:
            child_index += 1
        gc = self.forward_gc_map(obj, forced_index)
        return (
            tuple(None for _ in range(child_index))
            + (children[child_index].forward_map.map_gc(gc),)
            + tuple(None for _ in range(child_index, len(children) - 1))
        )

    def __str__(self) -> str:
        return "requirement placement strategy"

    def __repr__(self) -> str:
        return (
            f"RequirementPlacementStrategy(gcs={self.gcs}, "
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
        d["gcs"] = tuple(gc.to_jsonable() for gc in self.gcs)
        d["indices"] = self.indices
        d["direction"] = self.direction
        d["own_col"] = self.own_col
        d["own_row"] = self.own_row
        d["include_empty"] = self.include_empty
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "RequirementPlacementStrategy":
        gcs = tuple(GriddedChord.from_dict(gc) for gc in d.pop("gcs"))
        return cls(gcs=gcs, **d)
    
   
# Should this be renamed? It currently is "requirement placement" but perhaps it should be "chord placement"
class RequirementPlacementFactory(StrategyFactory[Tiling]):
    def __init__(self, max_reqlist_size: int = 1, 
                 max_chord_size: int = 1, 
                 max_width: int = 6, 
                 max_len: int = 6, 
                 dirs: Iterable[int] = (DIR_SOUTH, DIR_NORTH)) -> None:
        self.max_reqlist_size = max_reqlist_size
        self.max_chord_size = max_chord_size
        self.max_width = max_width
        self.max_len = max_len
        self.dirs = dirs

    def __call__(self, comb_class: Tiling) -> Iterator[RequirementPlacementStrategy]:
        reqs = comb_class.requirements
        for req_list in reqs:
            if len(req_list) <= self.max_reqlist_size and (
                all(len(req) <= self.max_chord_size for req in req_list)):
                for dir in self.dirs:
                    yield self._build_strategy(req_list, dir)

    def _build_strategy(self, req_list_to_place, dir):
        return RequirementPlacementStrategy(req_list_to_place, dir)
    
    def __str__(self) -> str:
        # sTODO fix this string
        s = f"Placing requirments lists up to size {self.max_reqlist_size}"
        s += f"with chords up to size {self.max_chord_size}"
        return s

    def __repr__(self) -> str:
        args = ", ".join(
            [
                f"max_reqlist_size={self.max_reqlist_size}",
                f"max_chord_size={self.max_chord_size}"
            ]
        )
        return f"{self.__class__.__name__}({args})"

    def to_jsonable(self) -> dict:
        d: dict = super().to_jsonable()
        interleaving = None

        d["max_reqlist_size"] = self.max_reqlist_size
        d["max_chord_size"] = self.max_chord_size
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "RequirementPlacementFactory":
        return cls(**d)

