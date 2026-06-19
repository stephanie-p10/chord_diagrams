import sys
from pathlib import Path
import time

# Allow running this example directly from the repo without requiring an
# editable install. (Adds `<repo>/src` to the import path.)
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from comb_spec_searcher import (
    AtomStrategy,
    CartesianProductStrategy,
    CombinatorialClass,
    CombinatorialObject,
    CombinatorialSpecificationSearcher,
    DisjointUnionStrategy,
    StrategyPack,
)
from src.common.chords import Chord, GriddedChord
from src.common import DIR_EAST, DIR_SOUTH
from src.strategies.chord_placement import ChordPlacementFactory
from src.strategies.factor import FactorFactory
from src.strategies.obstruction_inferral import (
    ObstructionInferralFactory,
    SubobstructionInferralFactory,
)
from src.strategies.requirement_insertion import RequirementInsertionFactory
from src.strategies.row_col_separation import RowColumnSeparationStrategy
from src.common.tiling import Tiling
start_time = time.time()



pack = StrategyPack(
    initial_strats=[FactorFactory(), ], # rules for initial strategies that are domain specific to be applied right away (like factor strategy)
    inferral_strats=[RowColumnSeparationStrategy(), SubobstructionInferralFactory()], # rules equivalence strategies (that get applied right away)  
    expansion_strats=[[RequirementInsertionFactory(maxreqlen=1, limited_insertion=1), ], [ChordPlacementFactory(max_reqlist_size=2, max_gc_size=3, dirs=(DIR_SOUTH, DIR_EAST))]], # rules for domain specific strategies that are used less often
    ver_strats=[AtomStrategy()], # returns a rule when the count of a class is known.
    name=("Finding specification for pattern avoiding chord diagrams (ex. non crossing)"),
)
# AtomStrategy is included in comb_spec_searcher, requires that is_atom and 
# minimum_size_of_object have been implemented in class of interest

non_crossing = Tiling(obstructions=(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),))

atom = Tiling(obstructions=(GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),
                            GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (0, 0), (0, 0))),
                            GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))),
              requirements=((GriddedChord(Chord((0, 0)), ((0,0), (0,0))),),))


empty = Tiling(obstructions=(GriddedChord(Chord((0,0)), ((0,0), (0,0))),))

contradiction = Tiling(obstructions=(GriddedChord(Chord((0, 0)), ((0,0), (0,0))),),
                       requirements=((GriddedChord(Chord((0,0)), ((0,0), (0,0))),),))

av012102 = Tiling(obstructions=(GriddedChord(Chord((0, 1, 2, 1, 0, 2)), ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),))


theorem303 = Tiling(obstructions=(GriddedChord(Chord((0, 1, 2, 0, 1, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
                                  GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
                                  GriddedChord(Chord((0, 1, 2, 1, 0, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),))

fig3_6 = Tiling(obstructions=(GriddedChord(Chord((0, 1, 2, 0, 1, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
                              GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),))

non_nesting = Tiling(obstructions=(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0),)*4),))

face_ex = Tiling((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),
                  GriddedChord(Chord((0, 1, 0, 1)), ((1, 0), (1, 0), (1, 0), (1, 0))),
                  GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))))

start_time = time.time()
searcher = CombinatorialSpecificationSearcher(fig3_6, pack)
spec = searcher.auto_search()

end_time = time.time()

gf = spec.get_genf()
print(gf)
print(end_time - start_time)

try:
    spec.show()
except Exception as e:
    # Some environments (headless CI, missing browser, restricted automation)
    # can't open the rendered HTML. The example still succeeds without it.
    print(f"spec.show() failed: {e!r}")


