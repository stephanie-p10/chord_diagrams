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
from strategies.chord_placement import RequirementPlacementStrategy

pack = StrategyPack(
    initial_strats=[FactorFactory(), RequirementInsertionFactory(maxreqlen=4)], # rules for initial strategies that are domain specific to be applied right away (like factor strategy)
    inferral_strats=[RowColumnSeparationStrategy(), RequirementPlacementStrategy([GriddedChord(Chord((0, 0)), ((0, 0), (0, 0)))], [0], DIR_SOUTH), ObstructionInferralFactory()], # rules equivalence strategies (that get applied right away)
    expansion_strats=[[]], # rules for domain specific strategies that are used less often
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

searcher = CombinatorialSpecificationSearcher(theorem303, pack)

spec = searcher.auto_search()
spec.get_genf()
spec.show()


#print(empty.is_atom())
#print(atom.is_atom())
#print(contradiction.is_atom())
#print(non_crossing.is_atom())



#print()
#print(empty.is_empty())
#print(atom.is_empty())
#print(contradiction.is_empty())
