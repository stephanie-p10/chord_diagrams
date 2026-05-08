import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

import time
from comb_spec_searcher import (
    AtomStrategy,
    CombinatorialSpecificationSearcher,
    StrategyPack,
)
from misc import DIR_EAST, DIR_SOUTH

from chords import GriddedChord, Chord
from tiling import Tiling
from strategies.obstruction_inferral import SubobstructionInferralFactory
from strategies.factor import FactorFactory
from strategies.requirement_insertion import RequirementInsertionFactory
from strategies.row_col_separation import RowColumnSeparationStrategy
from strategies.chord_placement import RequirementPlacementFactory


pack = StrategyPack(
    initial_strats=[FactorFactory(), ], # rules for initial strategies that are domain specific to be applied right away (like factor strategy)
    inferral_strats=[RowColumnSeparationStrategy(), SubobstructionInferralFactory()], # rules equivalence strategies (that get applied right away)  
    expansion_strats=[[RequirementInsertionFactory(maxreqlen=1, limited_insertion=1)], [RequirementPlacementFactory(max_reqlist_size=2, max_chord_size=3, dirs=(DIR_SOUTH, DIR_EAST))]], # rules for domain specific strategies that are used less often
    ver_strats=[AtomStrategy()], # returns a rule when the count of a class is known.
    name=("Finding specification for pattern avoiding chord diagrams (ex. non crossing)"),
)

non_crossing = Tiling(obstructions=(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),))

searcher = CombinatorialSpecificationSearcher(non_crossing, pack)
spec = searcher.auto_search()

print(spec.get_genf())
print(spec.get_terms(10))
spec.show()
