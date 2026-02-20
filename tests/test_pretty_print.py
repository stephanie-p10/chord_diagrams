import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from chords import GriddedChord, Chord
from tiling import Tiling

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