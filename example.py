import time
from comb_spec_searcher import (
    AtomStrategy,
    CartesianProductStrategy,
    CombinatorialClass,
    CombinatorialObject,
    CombinatorialSpecificationSearcher,
    DisjointUnionStrategy,
    StrategyPack,
)
from misc import DIR_EAST, DIR_SOUTH

from chords import GriddedChord, Chord
from tiling import Tiling
from strategies.obstruction_inferral import ObstructionInferralFactory
from strategies.factor import FactorFactory
from strategies.requirement_insertion import RequirementInsertionFactory
from strategies.row_col_separation import RowColumnSeparationStrategy
from strategies.chord_placement import RequirementPlacementFactory
start_time = time.time()


pack = StrategyPack(
    initial_strats=[FactorFactory(), RequirementInsertionFactory(maxreqlen=1, limited_insertion=False)], # rules for initial strategies that are domain specific to be applied right away (like factor strategy)
    inferral_strats=[RowColumnSeparationStrategy(), ObstructionInferralFactory()], # rules equivalence strategies (that get applied right away)
    expansion_strats=[[RequirementPlacementFactory(max_reqlist_size=2, max_chord_size=3, dirs=(DIR_SOUTH, DIR_EAST))]], # rules for domain specific strategies that are used less often
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


theorem303 = Tiling(obstructions=(GriddedChord(Chord((0, 1, 2, 0, 1, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
                                  GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
                                  GriddedChord(Chord((0, 1, 2, 1, 0, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),))

fig3_6 = Tiling(obstructions=(GriddedChord(Chord((0, 1, 2, 0, 1, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
                              GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),))

'''size_five_chords = Chord.of_length(5)
obs = [GriddedChord(chord, ((0, 0),)*10) for chord in size_five_chords]

size_four_atom = Tiling(obstructions=obs,
                       requirements=((GriddedChord(Chord((0, 1, 2, 3, 0, 1, 2, 3)), ((0, 0),)*8),),))'''

non_nesting = Tiling(obstructions=(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0),)*4),))

start_time = time.time()
searcher = CombinatorialSpecificationSearcher(fig3_6, pack)
spec = searcher.auto_search()

end_time = time.time()

gf = spec.get_genf()
print(gf)
print(end_time - start_time)

#spec.show()


