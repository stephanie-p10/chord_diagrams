import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from chords import GriddedChord, Chord
from tiling import Tiling
from algorithms.factor import Factor
from strategies.factor import FactorFactory, FactorStrategy
from strategies.requirement_insertion import RequirementInsertionStrategy
from strategies.chord_placement import RequirementPlacementStrategy



def crossed_chord(pos):
    return GriddedChord(Chord((0, 1, 0, 1)), (pos,) *4)
def single_chord(pos_left, pos_right):
    return GriddedChord(Chord((0,0)), (pos_left, pos_right))

non_crossing = Tiling(obstructions=(single_chord((1, 1), (3, 1)), 
                                    crossed_chord((1, 1)), 
                                    crossed_chord((3, 1)), 
                                    single_chord((0, 0), (0, 0)), 
                                    single_chord((2, 0), (2, 0)),
                                    GriddedChord(Chord((0,1)), ((0,0), (0,0))), 
                                    GriddedChord(Chord((1,0)), ((0,0), (0,0))),
                                    GriddedChord(Chord((0,1)), ((2,0), (2,0))), 
                                    GriddedChord(Chord((1,0)), ((2,0), (2,0)))),
                      requirements=([single_chord((0, 0), (2, 0))],))
#non_crossing.pretty_print_latex("non_crossing.tex")


# factory = FactorFactory(True)
# for strat in factory(non_crossing):
#     print("next strat:")
#     counter=1
#     for part in strat.decomposition_function(non_crossing):
#         part.pretty_print_latex("part"+str(counter)+".tex")
#         print("LEVEL 1, part "+str(counter)+":")
#         print(part)
#         for strat in factory(part):
#             print("Again:")
#             counter2=1
#             for part2 in strat.decomposition_function(part):
#                 print("LEVEL 2, part "+str(counter)+"_"+str(counter2)+":")
#                 print(part2)
#                 part2.pretty_print_latex("part"+str(counter)+"_"+str(counter2)+".tex")
#                 counter2+=1
#         counter+=1
#     print()


# first pattern in theorem 3.1.2 class
c1 = Chord((0, 1, 2, 0, 1, 2))
# second pattern in theorem 3.1.2 class
c2 = Chord((0, 1, 2, 0, 2, 1))
single_chord = Chord((0, 0))
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

t3_prime_to_t4_t5 = RequirementInsertionStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),))
t5_to_t5_prime = RequirementPlacementStrategy((GriddedChord(single_chord, ((1, 1), (3, 1))),), 0)
lukas_t5 = t3_prime_to_t4_t5.decomposition_function(lukas_t3_prime)[1]
lukas_t5._simplify()

lukas_t5_prime = t5_to_t5_prime.decomposition_function(lukas_t5)[0]

lukas_t5_prime.pretty_print_latex("lukas_t5_prime.tex")

factory = FactorFactory(True)
for strat in factory(lukas_t5_prime):
    print("next strat:")
    counter=1
    for part in strat.decomposition_function(lukas_t5_prime):
        part.pretty_print_latex("part"+str(counter)+".tex")
        print("LEVEL 1, part "+str(counter)+":")
        print(part)
        counter+=1
    print()