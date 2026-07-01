import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

import time
from comb_spec_searcher import (
    AtomStrategy,
    CombinatorialSpecificationSearcher,
    StrategyPack,
)
from src.common import DIR_EAST, DIR_SOUTH

from src.common.chords import GriddedChord, Chord
from src.common.tiling import Tiling
from src.strategies.obstruction_inferral import SubobstructionInferralFactory
from src.strategies.factor import FactorFactory
from src.strategies.requirement_insertion import RequirementInsertionFactory
from src.strategies.row_col_separation import RowColumnSeparationStrategy
from src.strategies.chord_placement import ChordPlacementFactory


pack = StrategyPack(
    initial_strats=[FactorFactory(), ], # rules for initial strategies that are domain specific to be applied right away (like factor strategy)
    inferral_strats=[RowColumnSeparationStrategy(), SubobstructionInferralFactory()], # rules equivalence strategies (that get applied right away)  
    expansion_strats=[[RequirementInsertionFactory(maxreqlen=1, limited_insertion=1)], [ChordPlacementFactory(max_reqlist_size=2, max_gc_size=3, dirs=(DIR_SOUTH, DIR_EAST))]], # rules for domain specific strategies that are used less often
    ver_strats=[AtomStrategy()], # returns a rule when the count of a class is known.
    name=("Finding specification for pattern avoiding chord diagrams (ex. non crossing)"),
)

non_crossing = Tiling(obstructions=(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),))

searcher = CombinatorialSpecificationSearcher(non_crossing, pack)
spec = searcher.auto_search()

print(spec.get_genf())
print(spec.get_terms(10))
#spec.show()

