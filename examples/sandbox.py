from cProfile import run
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.common.chords import Chord, GriddedChord
from src.common.tiling import Tiling
from itertools import product

import time

'''gc = GriddedChord(Chord((0, 0)), ((3, 1), (3, 1)))
print(gc.get_bounding_indices(1, 1))
'''
c1 = Chord((0, 1, 2, 0, 1, 2))
c2 = Chord((0, 1, 2, 0, 2, 1))
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

t_plain = Tiling(requirements=((GriddedChord(Chord((0,0)), ((0, 0), (0, 0))),),))

'''time1_start = time.time()
old_method = t5.all_chords_on_tiling(6)
time1_end = time.time()
print("old:")
for chord in old_method:
    print(chord)



time2_start = time.time()
new_method = t5.build_all_chords_on_tiling(3)
time2_end = time.time()
print("new")
for chord in new_method:
    print(chord)

print(time1_end-time1_start, time2_end-time2_start)


assert set(old_method) == set(new_method)

chords_built = t_plain.build_all_chords_on_tiling(3)

for chord in chords_built:
    print(chord)'''

print(t5)




