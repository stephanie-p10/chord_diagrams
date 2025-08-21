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
# sToDo: need factory for picking gcs to feed to strategy
class RequirementInsertionStrategy(DisjointUnionStrategy[Tiling, GriddedChord]):
    def __init__(self, gcs: Iterable[GriddedChord], ignore_parent: bool = False):
        super().__init__(ignore_parent=ignore_parent)
        self.gcs = frozenset(gcs)

    def decomposition_function(self, comb_class: Tiling) -> Tuple[Tiling, Tiling]:
        """
        Return a tuple of tiling. The first one avoids all the pattern in the
        list while the other contain one of the patterns in the list.
        """
        #print(comb_class.add_obstructions(self.gcs).obstructions)
        return comb_class.add_obstructions(self.gcs), comb_class.add_list_requirement(self.gcs)
    
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
    
    def backward_map(
        self,
        comb_class: Tiling,
        objs: Tuple[Optional[GriddedChord], ...],
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Iterator[GriddedChord]:
        if children is None:
            children = self.decomposition_function(comb_class)
        idx = DisjointUnionStrategy.backward_map_index(objs)
        return children[idx].backward_map.map_gc(cast(GriddedChord, objs[idx]))

    def forward_map(
        self,
        comb_class: Tiling,
        obj: GriddedChord,
        children: Optional[Tuple[Tiling, ...]] = None,
    ) -> Tuple[Optional[GriddedChord], Optional[GriddedChord]]:
        if children is None:
            children = self.decomposition_function(comb_class)
        if obj.avoids(*self.gcs):
            #print(children[0].forward_map)
            #print(children[0].forward_map.map_gc(obj))
            return (children[0].forward_map.map_gc(obj), None)
        return (None, children[1].forward_map.map_gc(obj))
    
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
        gcs = [GriddedChord.from_dict(gc) for gc in d.pop("gps")]
        return cls(gcs=gcs, **d)
    

class AbstractRequirementInsertionFactory(StrategyFactory[Tiling]):
    """
    Bases class for requirement insertion on tilings.

    It will create batch rules based on the containment or not of
    requirements.  The requirement  used are  provided by `req_list_to_insert`
    """

    def __init__(self, ignore_parent: bool = False):
        self.ignore_parent = ignore_parent

    @abc.abstractmethod
    def req_lists_to_insert(self, tiling: Tiling) -> Iterator[ListRequirement]:
        """
        Iterator over all the requirement list to insert to create the batch
        rules.
        """

    def __call__(self, comb_class: Tiling) -> Iterator[RequirementInsertionStrategy]:
        """
        Iterator over all the requirement insertion rules.
        """
        for req_list in list(set(self.req_lists_to_insert(comb_class))):
            yield RequirementInsertionStrategy(req_list, self.ignore_parent)

    def to_jsonable(self) -> dict:
        d: dict = super().to_jsonable()
        d["ignore_parent"] = self.ignore_parent
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "AbstractRequirementInsertionFactory":
        return cls(**d)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(ignore_parent={self.ignore_parent})"

class RequirementInsertionWithRestrictionFactory(AbstractRequirementInsertionFactory):
    """
    As RequirementInsertion, but a set of pattern to avoids and a maximum
    length can be provided.
    """

    def __init__(
        self,
        maxreqlen: int,
        extra_basis: Optional[List[Chord]] = None,
        ignore_parent: bool = False,
    ):
        if extra_basis is None:
            self.extra_basis = []
        else:
            assert isinstance(extra_basis, list) # EXTRA_BASIS_ERR
            assert all(isinstance(chord, Chord) for chord in extra_basis) # EXTRA_BASIS_ERR
            self.extra_basis = extra_basis
        self.maxreqlen = maxreqlen
        super().__init__(ignore_parent)

    def to_jsonable(self) -> dict:
        d: dict = super().to_jsonable()
        d["maxreqlen"] = self.maxreqlen
        d["extra_basis"] = self.extra_basis
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "RequirementInsertionWithRestrictionFactory":
        if d["extra_basis"] is None:
            extra_basis = None
        else:
            extra_basis = [Chord(p) for p in d["extra_basis"]]
        d.pop("extra_basis")
        return cls(extra_basis=extra_basis, **d)

    def __repr__(self) -> str:
        args = ", ".join(
            [
                f"maxreqlen={self.maxreqlen}",
                f"extra_basis={self.extra_basis}",
                f"ignore_parent={self.ignore_parent}",
            ]
        )
        return f"{self.__class__.__name__}({args})"
    
# sToDo: do we want all patts or just chord patts inserted?
class RequirementInsertionFactory(RequirementInsertionWithRestrictionFactory):
    """
    Insert all possible requirements the obstruction allows if the tiling does
    not have requirements.

    If <limited_insertion> is true, the default behavior, requirements will only be
    inserted on Tilings that have no requirements.
    """

    def __init__(
        self,
        maxreqlen: int = 2,
        extra_basis: Optional[List[Chord]] = None,
        limited_insertion: bool = True,
        ignore_parent: bool = False,
        allow_factorable_insertions: bool = False,
    ) -> None:
        self.limited_insertion = limited_insertion
        self.allow_factorable_insertions = allow_factorable_insertions
        super().__init__(maxreqlen, extra_basis, ignore_parent)

    def req_lists_to_insert(self, tiling: Tiling) -> Iterator[ListRequirement]:
        obs_tiling = Tiling(
            tiling.obstructions,
            remove_empty_rows_and_cols=False,
            derive_empty=False,
            simplify=False,
            #sorted_input=True,
        )
        """Generates lists to insert
        
        previously, this worked by finding natural number patterns up to length number of points.
        However, once these patterns are expanded with the Expansion strategy, it basically gives
        what would have been generated with chord diagrams
        
        A future TODO might be to figure out if it is better to add requirements based off of non 
        chord patterns or if chord patterns only should be used (which is the current behaviour)
        For example, non chord patterns would allow 01 = {0110, 0101, 0011} to be added as a reqlist, 
        but right now it only adds sigleton chord diagram requirement sets to be added."""
        for length in range(1, self.maxreqlen + 1):
            for gc in obs_tiling.all_chords_on_tiling(2 * length):
                if (self.allow_factorable_insertions or len(gc.factors()) == 1) and all(
                    p not in gc.patt for p in self.extra_basis
                ):
                    yield (GriddedChord(Chord(gc.patt), gc.pos),)

    def __call__(self, comb_class: Tiling) -> Iterator[RequirementInsertionStrategy]:
        if self.limited_insertion and comb_class.requirements:
            return
        yield from super().__call__(comb_class)

    def __str__(self) -> str:
        if self.maxreqlen == 1:
            return "point insertion"
        if self.extra_basis:
            perm_class = f" from Av{(self.extra_basis)}"
        else:
            perm_class = ""
        return f"requirement insertion{perm_class} up to length {self.maxreqlen}"

    def __repr__(self) -> str:
        args = ", ".join(
            [
                f"maxreqlen={self.maxreqlen}",
                f"extra_basis={self.extra_basis!r}",
                f"limited_insertion={self.limited_insertion}",
                f"ignore_parent={self.ignore_parent}",
                f"allow_factorable_insertions={self.allow_factorable_insertions}",
            ]
        )
        return f"{self.__class__.__name__}({args})"

    def to_jsonable(self) -> dict:
        d: dict = super().to_jsonable()
        d["limited_insertion"] = self.limited_insertion
        d["allow_factorable_insertions"] = self.allow_factorable_insertions
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "RequirementInsertionWithRestrictionFactory":
        if d["limited_insertion"] is None:
            extra_basis = None
        else:
            extra_basis = [Chord(c) for c in d["extra_basis"]]
        d.pop("extra_basis")
        limited_insertion = d.pop("limited_insertion")
        maxreqlen = d.pop("maxreqlen")
        allow_factorable_insertions = d.pop("allow_factorable_insertions", False)
        return cls(
            maxreqlen=maxreqlen,
            extra_basis=extra_basis,
            limited_insertion=limited_insertion,
            allow_factorable_insertions=allow_factorable_insertions,
            **d,
        )
    
#sToDo: other reqins factories could be moved over from permutations tilings
