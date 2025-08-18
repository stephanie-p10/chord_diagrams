import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from algorithms.simplify import SimplifyObstructionsAndRequirements

from chords import GriddedChord, Chord
from tiling import Tiling

ob_containing_ob = SimplifyObstructionsAndRequirements((GriddedChord(Chord((0, 1, 0, 2, 1, 2)), ((0,0),)*6), GriddedChord(Chord((0, 1, 0, 1)), ((0,0),)*4)),
                                                   (), 
                                                   (1, 1))
ob_containing_ob.remove_redundant_obstructions()
assert ob_containing_ob.obstructions == (GriddedChord(Chord((0, 1, 0, 1)), ((0,0),)*4),)

req_containing_ob = SimplifyObstructionsAndRequirements((GriddedChord(Chord((0, 0)), ((0,0), (0,0))),),
                                            ((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),),), 
                                            (1,1))
req_containing_ob.simplify()
assert req_containing_ob.obstructions == (GriddedChord(Chord((0, 0)), ((0,0), (0,0))),)
assert req_containing_ob.requirements == ((),)

"""all_from_21 = SimplifyObstructionsAndRequirements((GriddedChord(Chord((0,)), ((1, 0),)),
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

all_from_21.simplify()"""

req_containing_req = SimplifyObstructionsAndRequirements((),
                                                    ((GriddedChord(Chord((0, 1, 0, 2, 1, 2)), ((0, 0), )*6),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), )*4),),), 
                                                    (1,1))
req_containing_req.remove_redundant_requirements()
assert req_containing_req.requirements == ((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),),)

reqlist_containing_reqlist = SimplifyObstructionsAndRequirements((),
                                                    ((GriddedChord(Chord((0, 1, 0, 2, 1, 2)), ((0, 0), )*6),),
                                                      (GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), )*4),),), 
                                                    (1,1))
reqlist_containing_reqlist.remove_redundant_lists_requirements()
assert reqlist_containing_reqlist.requirements == (((GriddedChord(Chord((0, 1, 0, 2, 1, 2)), ((0, 0), )*6),),))

reqlist_containing_ob = SimplifyObstructionsAndRequirements((GriddedChord(Chord((0, 1, 0, 1)), ((0,0),)*4),),
                                                            ((GriddedChord(Chord((0, 1, 0, 1)), ((0,0),)*4),
                                                             GriddedChord(Chord((0, 1, 1, 0)), ((0,0),)*4),),),
                                                            (1, 1))
reqlist_containing_ob.remove_redundant_requirements()
assert reqlist_containing_ob.requirements == ((GriddedChord(Chord((0,1,1,0)), ((0,0),)*4),),)
assert reqlist_containing_ob.obstructions == (GriddedChord(Chord((0, 1, 0, 1)), ((0,0),)*4),)
