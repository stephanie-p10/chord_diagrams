import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

import abc
from itertools import chain, product
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple, cast

#import tilings.strategies as strat
from comb_spec_searcher import DisjointUnionStrategy, StrategyFactory
from comb_spec_searcher.exception import StrategyDoesNotApply
from comb_spec_searcher.strategies import Rule
from comb_spec_searcher.strategies.strategy import VerificationStrategy
from chords import Chord, GriddedChord
from tiling import Tiling
from tilings.algorithms import Factor, SubobstructionInferral

Cell = Tuple[int, int]
ListRequirement = Tuple[GriddedChord, ...]

"""
Requirement insertion splits a tiling into two disjoint sets.
For some set H of grids, a tiling T gets split into:
    Tiling T_1 which avoids all patterns in H
    Tiling T_2 which contains H
"""

# gcs = H
class RequirementInsertionStrategy(DisjointUnionStrategy[Tiling, GriddedChord]):
    def __init__(self, gcs: Iterable[GriddedChord], ignore_parent: bool = False):
        super().__init__(ignore_parent=ignore_parent)
        self.gcs = frozenset(gcs)

    def decomposition_function(self, comb_class: Tiling) -> Tuple[Tiling, Tiling]:
        """
        Return a tuple of tiling. The first one avoids all the pattern in the
        list while the other contain one of the patterns in the list.
        """
        return comb_class.add_obstructions(self.gcs), comb_class.add_list_requirement(
            self.gcs
        )
    
    def formal_step(self) -> str:
        """
        Return the formal step for the insertion according to the req_list
        inserted.

        This needs to be redefined if you want to insert list requirement with
        more than one requirement.
        """
        if len(self.gcs) == 1:
            req = tuple(self.gcs)[0]
            if req.is_localized():
                return f"insert {req.patt} in cell {req.pos[0]}"
            return f"insert {req}"
        raise NotImplementedError
    
    #sTodo: what is the purpose of the following two methods??
    def forward_map(self, comb_class, obj, children = None):
        return super().forward_map(comb_class, obj, children)
    
    #???
    def backward_map(
        self,
        comb_class: Tiling,
        objs: Tuple[Optional[GriddedChord], ...],
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Iterator[GriddedChord]:
        
        """Maps from a child class from the strategy to... """
        if children is None:
            children = self.decomposition_function(comb_class)
        idx = DisjointUnionStrategy.backward_map_index(objs)
        yield children[idx].backward_map.map_gc(cast(GriddedChord, objs[idx]))

    
    def __repr__(self) -> str:
        args = ", ".join([f"gcs={self.gcs}", f"ignore_parent={self.ignore_parent}"])
        return f"{self.__class__.__name__}({args})"

    def __str__(self) -> str:
        return "requirement insertion"

    def to_jsonable(self) -> dict:
        """Return a dictionary form of the strategy."""
        d: dict = super().to_jsonable()
        d.pop("inferrable")
        d.pop("possibly_empty")
        d.pop("workable")
        d["gcs"] = [gc.to_jsonable() for gc in self.gcs]
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "RequirementInsertionStrategy":
        gps = [GriddedChord.from_dict(gp) for gp in d.pop("gps")]
        return cls(gps=gps, **d)