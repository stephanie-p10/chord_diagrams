import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from chords import GriddedChord, Chord
from tiling import Tiling
from algorithms.factor import Factor
from strategies.factor import FactorFactory, FactorStrategy

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
non_crossing.pretty_print_latex("non_crossing.tex")


factory = FactorFactory(True)
for strat in factory(non_crossing):
    print("next strat:")
    counter=1
    for part in strat.decomposition_function(non_crossing):
        part.pretty_print_latex("part"+str(counter)+".tex")
        print("LEVEL 1, part "+str(counter)+":")
        print(part)
        for strat in factory(part):
            print("Again:")
            counter2=1
            for part2 in strat.decomposition_function(part):
                print("LEVEL 2, part "+str(counter)+"_"+str(counter2)+":")
                print(part2)
                part2.pretty_print_latex("part"+str(counter)+"_"+str(counter2)+".tex")
                counter2+=1
        counter+=1
    print()
