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
from strategies.factor import FactorFactory, FactorStrategy, Factor
from strategies.requirement_insertion import RequirementInsertionFactory, RequirementInsertionStrategy
from strategies.row_col_separation import RowColumnSeparationStrategy
from strategies.chord_placement import RequirementPlacementFactory, RequirementPlacementStrategy

# first pattern in theorem 3.1.2 class
c1 = Chord((0, 1, 2, 0, 1, 2))
# second pattern in theorem 3.1.2 class
c2 = Chord((0, 1, 2, 0, 2, 1))
single_chord = Chord((0, 0))

t1 = Tiling(obstructions=(GriddedChord(c1, ((0, 0),)*6), GriddedChord(c2, ((0, 0),)*6)),)
t2 = Tiling(obstructions=(GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),))
t3 = Tiling(obstructions=(GriddedChord(c1, ((0, 0),)*6), GriddedChord(c2, ((0, 0),)*6)),
            requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),),))
t3_prime = Tiling(obstructions=(GriddedChord.single_chord(((0, 0), (0, 0))),
                                GriddedChord.single_chord(((2, 0), (2, 0))), 
                                GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0))),

                                GriddedChord(c1, ((0, 0), (1, 1), (1, 1), (2, 0), (3, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),

                                GriddedChord(c2, ((0, 0), (1, 1), (1, 1), (2, 0), (3, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),),
                  requirements=((GriddedChord.single_chord(((0, 0), (2, 0))),),))
t4 = Tiling(obstructions=(GriddedChord(c1, ((1, 1),)*6), 
                          GriddedChord(c2, ((1, 1),)*6), 
                          GriddedChord(c1, ((3, 1),)*6), 
                          GriddedChord(c2, ((3, 1),)*6),
                          GriddedChord(single_chord, ((1, 1), (3, 1)))),
            requirements=((GriddedChord(single_chord, ((0, 0), (2, 0))),),))
t4_prime = Tiling(obstructions=(GriddedChord(c1, ((1, 1),)*6), 
                                GriddedChord(c2, ((1, 1),)*6), 
                                GriddedChord(c1, ((3, 2),)*6), 
                                GriddedChord(c2, ((3, 2),)*6),),
                  requirements=((GriddedChord(single_chord, ((0, 0), (2, 0))),),))

t6 = Tiling(obstructions=(GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                          GriddedChord(Chord((0, 0)), ((1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0))),),
            requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))),),))

t5 = Tiling(obstructions=(GriddedChord.single_chord(((0, 0), (0, 0))),
                          GriddedChord.single_chord(((2, 0), (2, 0))), 
                          GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                          GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                          GriddedChord(c1, ((0, 0), (1, 1), (1, 1), (2, 0), (3, 1), (3, 1))),
                          GriddedChord(c2, ((0, 0), (1, 1), (1, 1), (2, 0), (3, 1), (3, 1))),
                          GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                          GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                          GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                          GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                          GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (3, 1), (3, 1))),
                          GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (3, 1), (3, 1))),
                          GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (3, 1), (3, 1), (3, 1))),
                          GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (3, 1), (3, 1), (3, 1))),
                          GriddedChord(c1, ((1, 1), (1, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                          GriddedChord(c2, ((1, 1), (1, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                          GriddedChord(c1, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                          GriddedChord(c2, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                          GriddedChord(c1, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                          GriddedChord(c2, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),),
                  requirements=((GriddedChord.single_chord(((0, 0), (2, 0))),),
                                (GriddedChord.single_chord(((1, 1), (3, 1))),)))

lukas_t3_prime = Tiling(obstructions=(GriddedChord.single_chord(((0, 0), (0, 0))),
                                GriddedChord.single_chord(((2, 0), (2, 0))), 
                                GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0))),

                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),

                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                
                                GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (3, 1), (3, 1))),
                                GriddedChord(Chord((0, 1, 1, 0)), ((1, 1), (1, 1), (3, 1), (3, 1))),),
                  requirements=((GriddedChord.single_chord(((0, 0), (2, 0))),),))




req_ins_t1_to_t2_t3 = RequirementInsertionStrategy((GriddedChord(single_chord, ((0, 0), (0, 0))),))
req_pl_t3 = RequirementPlacementStrategy((GriddedChord(single_chord, ((0, 0), (0 ,0))),), 3)
t3_prime_to_t4_t5 = RequirementInsertionStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),))
t5_to_t5_prime = RequirementPlacementStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),), 0)


t3_prime_to_t4_t5 = RequirementInsertionStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),))

t5_to_t5_prime = RequirementPlacementStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),), 0)

t5_prime = t5_to_t5_prime.decomposition_function(t5)[0]

lukas_t5 = t3_prime_to_t4_t5.decomposition_function(lukas_t3_prime)[1]

lukas_t5_prime = t5_to_t5_prime.decomposition_function(lukas_t5)[0]
print(lukas_t5_prime)

#for child in t5_to_t5_prime.decomposition_function(t5):
#    print(child)


for comp in Factor(lukas_t5_prime).get_components():
    print(comp)

#assert set(req_ins_t1_to_t2_t3.decomposition_function(t1)) == set([t2, t3]) # strategy 1
#assert set(req_pl_t3.decomposition_function(t3)) == set([t3_prime]) # strategy 2

#print("asserts passed")








