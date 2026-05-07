import _direct_run_bootstrap  

from chord_diagrams.chords import GriddedChord, Chord
from chord_diagrams.tiling import Tiling

tiling2 = Tiling(
        obstructions=[
            GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), ) * 4),
            GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
            GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))
        ],
        requirements=[
            [GriddedChord(Chord((0, 0)), ((1, 1), (2, 1)))],
            [GriddedChord(Chord((0, 1)), ((0, 0), (1, 1)))],
        ],
    )

print("Initialized tiling2")