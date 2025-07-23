import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from algorithms.obstruction_inferral import *
from chords import GriddedChord, Chord
from tiling import Tiling

from strategies.obstruction_inferral import ObstructionInferralFactory

def crossed_chord(pos):
    return GriddedChord(Chord((0, 1, 0, 1)), (pos,) *4)
def single_chord(pos_left, pos_right):
    return GriddedChord(Chord((0,0)), (pos_left, pos_right))
#print(crossed_chord((0,0)))

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

tiling1 = Tiling((single_chord((0,0),(2,0)), 
             single_chord((1,0), (2,0)), 
             single_chord((2,0), (2,0)), 
             GriddedChord(Chord((0,1)), ((2,0), (2,0))), 
             GriddedChord(Chord((1,0)), ((0,1), (1,1))), 
             GriddedChord(Chord((1,0)), ((2,0), (2,0))),), 
             ((single_chord((0,0), (1,0)),), 
              (single_chord((1,1), (2,1)),)))
tiling2 = Tiling((GriddedChord(Chord((0,1,0,1)), ((1,1),) *4),), 
                 ((single_chord((0,0), (2,0)),), (GriddedChord(Chord((0,1)), ((0,0), (1,1))),)))
tilings = [non_crossing, tiling1, tiling2]

subobs_nc = SubobstructionInferral(non_crossing)
subobs_t1 = SubobstructionInferral(tiling2)
subobs_t2 = SubobstructionInferral(tiling2)

all_obs_nc = AllObstructionInferral(non_crossing, 2)
all_obs_t1 = AllObstructionInferral(tiling1, 4)
all_obs_t2 = AllObstructionInferral(tiling2, 4)

emptycell_nc = EmptyCellInferral(non_crossing)
emptycell_t1 = EmptyCellInferral(tiling1)
emptycell_t2 = EmptyCellInferral(tiling2)

#for chord in non_crossing.all_chords_on_tiling(4, True):
   # print(chord)

#print()

assert set(all_obs_nc.new_obs()) == set([GriddedChord(Chord((0, 1)), ((0, 0), (2, 0))),
                                         GriddedChord(Chord((1, 0)), ((0, 0), (2, 0))),
                                         GriddedChord(Chord((1, 0)), ((1, 1), (3, 1)))])

assert set(subobs_nc.new_obs()) == set([GriddedChord(Chord((1, 0, 1)), ((1, 1), (1, 1), (1, 1))),
                                         GriddedChord(Chord((1, 0, 1)), ((3, 1), (3, 1), (3, 1)))])


factory = ObstructionInferralFactory()

for strat in factory(non_crossing):
    for item in strat.decomposition_function(non_crossing):
        print(item)

