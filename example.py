from comb_spec_searcher import (
    AtomStrategy,
    CartesianProductStrategy,
    CombinatorialClass,
    CombinatorialObject,
    CombinatorialSpecificationSearcher,
    DisjointUnionStrategy,
    StrategyPack,
)

from strategies.requirement_insertion import ReqIns

from chords import GriddedChord, Chord
from tiling import Tiling

pack = StrategyPack(
    initial_strats=[], # rules for initial strategies that are domain specific to be applied right away (like factor strategy)
    inferral_strats=[], # rules equivalence strategies (that get applied right away)
    expansion_strats=[[]], # rules for domain specific strategies that are used less often
    ver_strats=[AtomStrategy()], # returns a rule when the count of a class is known.
    name=("Finding specification for pattern avoiding chord diagrams (ex. non crossing)"),
)

# AtomStrategy is included in comb_spec_searcher, requires that is_atom and 
# minimum_size_of_object have been implemented in class of interest



