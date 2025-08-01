import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from algorithms.simplify import SimplifyObstructionsAndRequirements

from chords import GriddedChord, Chord
from tiling import Tiling

obs_contains = SimplifyObstructionsAndRequirements((GriddedChord(Chord((0, 1, 0, 1)), ((0,0),)*4), GriddedChord(Chord((0, 0)), ((0,0), (0,0)))),
                                                   (), 
                                                   (1, 1))

obs_contains.remove_redundant_obstructions()

assert obs_contains.obstructions == (GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),)

size_2_obs = SimplifyObstructionsAndRequirements((GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                                                  GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                                                  GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                                                  ),
                                                  (), 
                                                  (2, 1))

size_2_obs.remove_redundant_obstructions()

#print(size_2_obs.obstructions)

empty = SimplifyObstructionsAndRequirements((GriddedChord(Chord((0, 0)), ((0,0), (0,0))),),
                                            ((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),),), 
                                            (1,1))
print(empty.obstructions, empty.requirements)

empty.simplify()

print(empty.obstructions, empty.requirements)

all_from_21 = SimplifyObstructionsAndRequirements((GriddedChord(Chord((0,)), ((1, 0),)),
                                                   GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                                                   GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))),

                                                   GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (2, 0))),
                                                   GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                                   GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (2, 0), (2, 0), (2, 0))),

                                                   GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (2, 0))),
                                                   GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                                   GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (2, 0), (2, 0), (2, 0))),

                                                   GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                                   GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (2, 0), (2, 0), (2, 0))),

                                                    GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))),
                                                    GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))),
                                                    GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))),
                                                    GriddedChord(Chord((1, 0)), ((2, 0), (2, 0)))),
                                                    ((GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),),),
                                                    (3, 1))

all_from_21.simplify()

print()
#print(all_from_21.obstructions)
print()
#print(all_from_21.requirements)
print()

redundant_req = SimplifyObstructionsAndRequirements((),
                                                    ((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),),
                                                     (GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),),), 
                                                    (1,1))

redundant_req.remove_redundant_lists_requirements()
print(redundant_req.requirements)
