import sys
from pathlib import Path
 
_src_root = Path(__file__).resolve().parents[2]  # .../src
if str(_src_root) not in sys.path:
    sys.path.insert(0, str(_src_root))

from comb_spec_searcher import (
    AtomStrategy,
    CartesianProductStrategy,
    CombinatorialClass,
    CombinatorialObject,
    CombinatorialSpecificationSearcher,
    DisjointUnionStrategy,
    StrategyPack,
)
from src.common import DIR_EAST, DIR_SOUTH

from src.common.chords import GriddedChord, Chord
from src.common.tiling import Tiling
from src.strategies.obstruction_inferral import SubobstructionInferralFactory
from src.strategies.factor import FactorFactory, FactorStrategy, Factor
from src.strategies.requirement_insertion import RequirementInsertionFactory, RequirementInsertionStrategy
from src.strategies.row_col_separation import RowColumnSeparationStrategy
from src.strategies.chord_placement import ChordPlacementFactory, RequirementPlacementStrategy, RequirementPlacement

# patterns from in theorem 3.1.2 class
c1 = Chord((0, 1, 2, 0, 1, 2))
c2 = Chord((0, 1, 2, 0, 2, 1))

t1 = Tiling(obstructions=(GriddedChord(c1, ((0, 0),)*6), GriddedChord(c2, ((0, 0),)*6)),)

t2 = Tiling(obstructions=(GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),))

t3 = Tiling(obstructions=(GriddedChord(c1, ((0, 0),)*6), GriddedChord(c2, ((0, 0),)*6)),
            requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),),))

strat_t1_to_t2_t3 = RequirementInsertionStrategy((GriddedChord.single_chord(((0, 0), (0, 0))),))

assert strat_t1_to_t2_t3.decomposition_function(t1) == (t2, t3)

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

strat_t3_to_t3p = RequirementPlacementStrategy((GriddedChord.single_chord(((0, 0), (0, 0))),), 3)

assert strat_t3_to_t3p.decomposition_function(t3) == (t3_prime,)

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

subobs_fact = SubobstructionInferralFactory()

strat_t3p_to_lukas_t3p = next(subobs_fact(t3_prime))

assert strat_t3p_to_lukas_t3p.decomposition_function(t3_prime) == (lukas_t3_prime,)

t4 = Tiling(obstructions=(GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                          GriddedChord(Chord((0, 0)), ((1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0)))),
            requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))),),))

t5 = Tiling(obstructions=(GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),
                          GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0))),
                          GriddedChord(c1, ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c1, ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),

                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),
                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0))),
                          GriddedChord(c2, ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c2, ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                                
                          GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0))),))

factor_fact = FactorFactory()

strat_lukas_t3p_to_t4_t5 = next(factor_fact(lukas_t3_prime))

assert strat_lukas_t3p_to_t4_t5.decomposition_function(lukas_t3_prime) == (t4, t5)

t6 = Tiling(obstructions=(GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),
                          GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0))),
                          GriddedChord(c1, ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c1, ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),

                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),
                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0))),
                          GriddedChord(c2, ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c2, ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                                
                          GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord.single_chord(((0, 0), (1, 0)))))

t7 = Tiling(obstructions=(GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),
                          GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0))),
                          GriddedChord(c1, ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c1, ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),

                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),
                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0))),
                          GriddedChord(c2, ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c2, ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                                
                          GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0)))),
            requirements=((GriddedChord.single_chord(((0, 0), (1, 0))),),))

strat_t5_to_t6_t7 = RequirementInsertionStrategy((GriddedChord.single_chord(((0, 0), (1, 0))),), 3)

assert strat_t5_to_t6_t7.decomposition_function(t5) == (t6, t7)

t6pp = Tiling(obstructions=(GriddedChord(c1, ((1, 0), )*6), 
                            GriddedChord(c2, ((1, 0), )*6),
                            GriddedChord(c1, ((0, 1), )*6), 
                            GriddedChord(c2, ((0, 1), )*6),))

strat_t6_to_t6pp = RowColumnSeparationStrategy()

assert strat_t6_to_t6pp.decomposition_function(t6) == (t6pp,)

strat_t6pp_to_t1_t1 = next(factor_fact(t6pp))

assert strat_t6pp_to_t1_t1.decomposition_function(t6pp) == (t1, t1)

strat_t7_to_t8 = RequirementPlacementStrategy([GriddedChord.single_chord(((0, 0), (1, 0)))], 3)

print(strat_t7_to_t8.decomposition_function(t7)[0])









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




e7 = Tiling(obstructions=(GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),
                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0))),
                          GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0))),
                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (1, 0))),
                          GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord(c1, ((0, 0), (0, 0), (0, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c2, ((0, 0), (0, 0), (0, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c1, ((0, 0), (0, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c2, ((0, 0), (0, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c1, ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c2, ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c1, ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),
                          GriddedChord(c2, ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))),),
            requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))),),))



'''req_pl_t3 = RequirementPlacementStrategy((GriddedChord(single_chord, ((0, 0), (0 ,0))),), 3)
t3_prime_to_t4_t5 = RequirementInsertionStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),))
t5_to_t5_prime = RequirementPlacementStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),), 0)


t3_prime_to_t4_t5 = RequirementInsertionStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),))

t5_to_t5_prime = RequirementPlacementStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),), 0)

t5_prime = t5_to_t5_prime.decomposition_function(t5)[0]



lukas_t5 = t3_prime_to_t4_t5.decomposition_function(lukas_t3_prime)[1]
lukas_t5._simplify()

req_pl_algo = RequirementPlacement(lukas_t5)

place_one_point = req_pl_algo.place_point(GriddedChord(Chord((0, 0)), ((1, 1), (3, 1))), 0, True)
place_one_point._simplify()

req_pl_point_two_algo = RequirementPlacement(place_one_point)
#print(place_one_point)

place_two_point = req_pl_point_two_algo.place_point(GriddedChord(Chord((0, 0)), ((2, 2), (5, 2))), 1, False)
place_two_point._simplify()

#print(place_two_point)

lukas_t5_prime = t5_to_t5_prime.decomposition_function(lukas_t5)[0]
#lukas_t5_prime._simplify()
#print(lukas_t5_prime)



#for child in t5_to_t5_prime.decomposition_function(t5):
#    print(child)''''''

non_crossing_ex1 = Tiling((GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (1, 1), (1, 1))),
                       GriddedChord(Chord((0, 1, 0, 1)), ((3, 1), (3, 1), (3, 1), (3, 1))), 
                       GriddedChord(Chord((0, 0)), ((1, 1), (3, 1))), 
                       GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))), 

                       GriddedChord(Chord((0, 1)), ((1, 1), (3, 1))), 
                       GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))), 
                       GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((1, 0)), ((2, 0), (2, 0))), ),
                       ((GriddedChord(Chord((0, 0)), ((0, 0), (3, 0))),),))

for comp in Factor(non_crossing_ex1).get_components():
    print(comp)

#assert set(req_ins_t1_to_t2_t3.decomposition_function(t1)) == set([t2, t3]) # strategy 1
#assert set(req_pl_t3.decomposition_function(t3)) == set([t3_prime]) # strategy 2

#print("asserts passed")

#print("t5_prime:")

#print(t5_prime)



# t1.pretty_print_latex("t1.tex")

# t2.pretty_print_latex("t2.tex")

t3.pretty_print_latex("t3.tex")

# t3_prime.pretty_print_latex("t3_prime.tex")

t4.pretty_print_latex("t4.tex")

# t4_prime.pretty_print_latex("t4_prime.tex")

t5.pretty_print_latex("t5.tex")

# t5_prime.pretty_print_latex("t5_prime.tex")

lukas_t5.pretty_print_latex("lukas_t5.tex")

lukas_t5_prime.pretty_print_latex("lukas_t5_prime.tex")'''