import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from chords import Chord, GriddedChord
from algorithms.expansion import Expansion

gc = GriddedChord(Chord((0, 1)), ((0, 0),)*2)

ex_01_algo = Expansion([], [], (1, 1))

gc_0_added = ex_01_algo._expand_single_point(gc, [])

assert ex_01_algo.expand_gridded_chord(gc_0_added[0], []) == [GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))), 
                                                                  GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))]
assert ex_01_algo.expand_gridded_chord(gc, []) == [GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (0, 0), (0, 0))), 
                                                       GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))), 
                                                       GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))]

isolated_chord_obs = [GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                      GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))), 
                      GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))), 
                      GriddedChord(Chord((0, 0)), ((1, 0), (1, 0))),
                      GriddedChord(Chord((0, 1)), ((1, 0), (1, 0))), 
                      GriddedChord(Chord((1, 0)), ((1, 0), (1, 0))),]

isolated_chord_algo = Expansion(isolated_chord_obs, ((GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))),),), (2, 1))
isolated_chord_algo.expand_obstructions()
assert set(isolated_chord_algo._obstructions) == set([GriddedChord(Chord((0, 0)), ((0 ,0), (0, 0))),
                                                      GriddedChord(Chord((0, 0)), ((1, 0), (1, 0))),
                                                      GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0)))])


point_obs = [GriddedChord(Chord((0,)), ((0, 0),)),
             GriddedChord(Chord((0,)), ((1, 0),)),
             GriddedChord(Chord((0,)), ((0, 1),)),
             GriddedChord(Chord((0,)), ((1, 1),))]
point_obs_algo = Expansion(point_obs, (), (2, 2))
point_obs_algo.expand_obstructions()
assert set(point_obs_algo._obstructions) == set([GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                                                 GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))),
                                                 GriddedChord(Chord((0, 0)), ((1, 0), (1, 0))),
                                                 GriddedChord(Chord((0, 0)), ((0, 1), (0, 1))),
                                                 GriddedChord(Chord((0, 0)), ((0, 1), (1, 1))),
                                                 GriddedChord(Chord((0, 0)), ((1, 1), (1, 1))),])


size_two_reqs = [GriddedChord(Chord((0, 1)), ((0, 0), (0, 0)))]
size_two_reqs_algo = Expansion((), (size_two_reqs,), (1, 1))
size_two_reqs_algo.expand_requirements()
assert set([tuple(list) for list in size_two_reqs_algo._requirements]) == set([(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (0, 0), (0, 0))), 
                                                     GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0) )), 
                                                     GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0) )),),])

print("asserts passed")


