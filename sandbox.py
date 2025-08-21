from chords import Chord, GriddedChord
from tiling import Tiling
from algorithms.expansion import Expansion

'''all_simplify = Tiling(obstructions= (GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                                     GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))),
                                     GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),
                                     GriddedChord(Chord((0, 0)), ((0, 0), (3, 0))),
                                     GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),
                                     GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                                     GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))),
                                     ),
                     requirements=(#(GriddedChord(Chord((0, 1, 2)), ((1, 1), (1, 1), (2, 2))),),
                                   (GriddedChord(Chord((0, 1)), ((1, 1), (1, 1))),), # is implied by req 3
                                   (GriddedChord(Chord((0, 1, 0, 1)), ((1, 1),) * 4),),
                                   (GriddedChord(Chord((0, 0)), ((1, 0), (2, 0))), GriddedChord(Chord((0, 1, 1, 0)), ((1, 0), (1, 0), (2, 0), (2, 0)))),
                                   ),)

print(all_simplify.requirements)
print(all_simplify.obstructions)'''
print(bool(Chord()))