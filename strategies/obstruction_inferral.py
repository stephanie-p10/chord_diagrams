import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, cast

from comb_spec_searcher import DisjointUnionStrategy, StrategyFactory
from comb_spec_searcher.exception import StrategyDoesNotApply

from chords import GriddedChord
from tiling import Tiling

class ObstructionInferralStrategy(DisjointUnionStrategy[Tiling, GriddedChord]):
    def __init__(self, gcs: Iterable[GriddedChord]):
        self.gcs = tuple(sorted(gcs))
        super().__init__(
            ignore_parent=True, inferrable=True, possibly_empty=False, workable=True
        )

    def decomposition_function(self, comb_class: Tiling) -> Tuple[Tiling]:
        return (comb_class.add_obstructions(self.gcs),)

    def formal_step(self) -> str:
        """Return a string describing the operation performed."""
        if all(gc.is_point() for gc in self.gcs):
            empty_cells_str = ", ".join(map(str, (gc.pos[0] for gc in self.gcs)))
            return f"the cells {{{empty_cells_str}}} are empty"
        return f"added the obstructions {{{', '.join(map(str, self.gcs))}}}"

    def extra_parameters(
        self, comb_class: Tiling, children: Optional[Tuple[Tiling, ...]] = None
    ) -> Tuple[Dict[str, str], ...]:
        if not comb_class.extra_parameters:
            return super().extra_parameters(comb_class, children)
        if children is None:
            children = self.decomposition_function(comb_class)
            if children is None:
                raise StrategyDoesNotApply("Strategy does not apply")
        child = children[0]
        params: Dict[str, str] = {}
        for assumption in comb_class.assumptions:
            mapped_assumption = (assumption).avoiding( # fix with child.forward_map.map_assumption when implemented
                child.obstructions
            )
            if mapped_assumption.gcs:
                parent_var = comb_class.get_assumption_parameter(assumption)
                child_var = child.get_assumption_parameter(mapped_assumption)
                params[parent_var] = child_var
        return (params,)

    def backward_map(
        self,
        comb_class: Tiling,
        objs: Tuple[Optional[GriddedChord], ...],
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Iterator[GriddedChord]:
        if children is None:
            children = self.decomposition_function(comb_class)
        yield children[0] # add .backward_map.map_gc(cast(GriddedChord, objs[0])) when implemented

    def forward_map(
        self,
        comb_class: Tiling,
        obj: GriddedChord,
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Tuple[GriddedChord]:
        if children is None:
            children = self.decomposition_function(comb_class)
        return (children[0],) # add .forward_map.map_gc(obj) when implemented

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(gcs={self.gcs})"

    def __str__(self) -> str:
        return self.formal_step()

    # JSON methods
    def to_jsonable(self) -> dict:
        """Return a dictionary form of the strategy."""
        d: dict = super().to_jsonable()
        d.pop("ignore_parent")
        d.pop("inferrable")
        d.pop("possibly_empty")
        d.pop("workable")
        d["gcs"] = [gc.to_jsonable() for gc in self.gcs]
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "ObstructionInferralStrategy":
        gcs = [GriddedChord.from_dict(gc) for gc in d.pop("gcs")]
        assert not d
        return cls(gcs=gcs)
