import _direct_run_bootstrap

from chord_diagrams.chords import GriddedChord, Chord
from chord_diagrams.tiling import Tiling

# tiling2 = Tiling(
#         obstructions=[
#             GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), ) * 4),
#             GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
#             GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))
#         ],
#         requirements=[
#         ],
#     )
atom = Tiling(obstructions=(GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),
                            GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (0, 0), (0, 0))),
                            GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))),
              requirements=((GriddedChord(Chord((0, 0)), ((0,0), (0,0))),),))

#atom.pretty_print_latex("atom.tex")

# non_crossing = Tiling(obstructions=(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),))

#non_crossing.pretty_print_latex("non_crossing.tex")

# theorem303 = Tiling(obstructions=(GriddedChord(Chord((0, 1, 2, 0, 1, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
#                                   GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),
#                                   GriddedChord(Chord((0, 1, 2, 1, 0, 2)), ((0, 0), (0,0), (0,0), (0,0), (0,0), (0,0))),))

# theorem303.pretty_print_latex("theorem303.tex")

place_point = Tiling(obstructions=(GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                                   GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))),
                                   GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                   GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 1), (2, 0), (3, 1))),
                                   GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (1, 1), (1, 1))),
                                   GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (1, 1), (3, 1))),
                                   GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (3, 1), (3, 1))),
                                   GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (3, 1), (3, 1), (3, 1))),
                                   GriddedChord(Chord((0, 1, 0, 1)), ((3, 1), (3, 1), (3, 1), (3, 1))),
                                   GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0)))),
                    requirements=((GriddedChord(Chord((0, 0)), ((0,0), (2,0))),),))

print(place_point)

place_point.pretty_print_latex("place_point.tex")

#print(atom.chord_cells)
print("Chord Row Cells:")
print(place_point.chord_row_cells)

atom.pretty_print_latex("atom.tex")

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