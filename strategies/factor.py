import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from chords import GriddedChord, Chord
from tiling import Tiling

from collections import Counter
from functools import reduce
from itertools import chain, combinations
from operator import mul
from typing import (
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    Set,
    Tuple,
    Union,
    cast,
)

from sympy import Eq, Function

from comb_spec_searcher import (
    CartesianProduct,
    CartesianProductStrategy,
    StrategyFactory,
)
from comb_spec_searcher.exception import StrategyDoesNotApply
from comb_spec_searcher.typing import (
    Parameters,
    SubObjects,
    SubRecs,
    SubSamplers,
    SubTerms,
    Terms,
)
from tilings.algorithms import (
    Factor,
    FactorWithInterleaving,
    FactorWithMonotoneInterleaving,
)
from tilings.assumptions import TrackingAssumption
from tilings.exception import InvalidOperationError
from tilings.misc import multinomial, partitions_iterator
from tilings.strategies.assumption_insertion import (
    AddAssumptionsConstructor,
    AddAssumptionsStrategy,
)

Cell = Tuple[int, int]

__all__ = (
    "FactorFactory",
    "FactorStrategy",
    "FactorWithInterleavingStrategy",
    "FactorWithMonotoneInterleavingStrategy",
)

TempGP = Tuple[
    Tuple[
        Tuple[Union[float, int], ...], Tuple[Union[float, int], ...], Tuple[Cell, ...]
    ],
    ...,
]


class FactorStrategy(CartesianProductStrategy[Tiling, GriddedChord]):
    def __init__(
        self,
        partition: Iterable[Iterable[Cell]],
        ignore_parent: bool = True,
        workable: bool = True,
    ):
        self.partition = tuple(sorted(tuple(sorted(p)) for p in partition))
        super().__init__(
            ignore_parent=ignore_parent, workable=workable, inferrable=False
        )

    def decomposition_function(self, comb_class: Tiling) -> Tuple[Tiling, ...]:
        return tuple(comb_class.sub_tiling(cells) for cells in self.partition)

    def extra_parameters(
        self, comb_class: Tiling, children: Optional[Tuple[Tiling, ...]] = None
    ) -> Tuple[Dict[str, str], ...]:
        if children is None:
            children = self.decomposition_function(comb_class)
            if children is None:
                raise StrategyDoesNotApply("Strategy does not apply")
        extra_parameters: Tuple[Dict[str, str], ...] = tuple({} for _ in children)
        for parent_var, assumption in zip(
            comb_class.extra_parameters, comb_class.assumptions
        ):
            for idx, child in enumerate(children):
                new_assumption = child.forward_map.map_assumption(assumption).avoiding(
                    child.obstructions
                )
                if new_assumption.gcs:
                    child_var = child.get_assumption_parameter(new_assumption)
                    extra_parameters[idx][parent_var] = child_var
        return extra_parameters

    def formal_step(self) -> str:
        """
        Return a string that describe the operation performed on the tiling.
        """
        partition_str = " / ".join(
            f"{{{', '.join(map(str, part))}}}" for part in self.partition
        )
        return f"factor with partition {partition_str}"

    def backward_map(
        self,
        comb_class: Tiling,
        objs: Tuple[Optional[GriddedChord], ...],
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Iterator[GriddedChord]:
        if children is None:
            children = self.decomposition_function(comb_class)
        gcs_to_combine = tuple(
            tiling.backward_map.map_gc(cast(GriddedChord, gc))
            for gc, tiling in zip(objs, children)
        )
        temp = [
            ((cell[0], idx), (cell[1], val))
            for gc in gcs_to_combine
            for (idx, val), cell in zip(enumerate(gc.patt), gc.pos)
        ]
        temp.sort()
        new_pos = [(idx[0], val[0]) for idx, val in temp]
        new_patt = Chord.to_standard(val for _, val in temp)
        yield GriddedChord(new_patt, new_pos)

    def forward_map(
        self,
        comb_class: Tiling,
        obj: GriddedChord,
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Tuple[GriddedChord, ...]:
        if children is None:
            children = self.decomposition_function(comb_class)
        return tuple(
            tiling.forward_map.map_gc(obj.get_subchord_in_cells(part)) # this only gives chords that have both endpoints in part
            for tiling, part in zip(children, self.partition)
        )

    def __str__(self) -> str:
        return "factor"

    def __repr__(self) -> str:
        args = ", ".join(
            [
                f"partition={self.partition}",
                f"ignore_parent={self.ignore_parent}",
                f"workable={self.workable}",
            ]
        )
        return f"{self.__class__.__name__}({args})"

    # JSON methods

    def to_jsonable(self) -> dict:
        """Return a dictionary form of the strategy."""
        d: dict = super().to_jsonable()
        d.pop("inferrable")
        d.pop("possibly_empty")
        d["partition"] = self.partition
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "FactorStrategy":
        partition = cast(
            Tuple[Tuple[Cell]],
            tuple(tuple(tuple(c) for c in p) for p in d.pop("partition")),
        )
        return cls(partition=partition, **d)
    


class FactorFactory(StrategyFactory[Tiling]):
    def __init__(
        self,
        unions: bool = False,
        ignore_parent: bool = True,
        workable: bool = True,
    ) -> None:
        
        self.factor_class = FactorStrategy
        self.unions = unions
        self.ignore_parent = ignore_parent
        self.workable = workable
        self.tracked = False

    def __call__(self, comb_class: Tiling) -> Iterator[FactorStrategy]:
        factor_algo = Factor(comb_class)
        print("factor_algo: ", factor_algo)
        if factor_algo.factorable():
            min_comp = tuple(tuple(part) for part in factor_algo.get_components())
            if self.unions:
                for partition in partitions_iterator(min_comp):
                    components = tuple(
                        tuple(chain.from_iterable(part)) for part in partition
                    )
                    yield self._build_strategy(components, workable=False)
            yield self._build_strategy(min_comp, workable=self.workable)

    def _build_strategy(
        self, components: Tuple[Tuple[Cell, ...], ...], workable: bool
    ) -> FactorStrategy:
        """
        Build the factor strategy for the given components.

        It ensure that a plain factor rule is returned.
        """
        return FactorStrategy(
            components, ignore_parent=self.ignore_parent, workable=workable
        )

    def __str__(self) -> str:
        s = "factor"
        if self.unions:
            s = "unions of " + s
        if self.tracked:
            s = "tracked " + s
        return s

    def __repr__(self) -> str:
        interleaving = None
        args = ", ".join(
            [
                f"unions={self.unions}",
                f"ignore_parent={self.ignore_parent}",
                f"workable={self.workable}",
            ]
        )
        return f"{self.__class__.__name__}({args})"

    def to_jsonable(self) -> dict:
        d: dict = super().to_jsonable()
        interleaving = None

        d["unions"] = self.unions
        d["ignore_parent"] = self.ignore_parent
        d["workable"] = self.workable
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "FactorFactory":
        return cls(**d)
