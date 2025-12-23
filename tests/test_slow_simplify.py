import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from chords import GriddedChord, Chord
from tiling import Tiling

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