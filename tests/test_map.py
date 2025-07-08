import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from chords import GriddedChord, Chord
from algorithms.map import RowColMap

id = RowColMap.identity((3,3))
double = RowColMap({0:0, 1:2, 2:4, 3:6}, {0:0, 1:2, 2:4}, False)
not_inj = RowColMap({0:0, 1:1, 2:2, 3:2}, {0:0, 1:1, 2:2, 3:3})

crossed = GriddedChord(Chord((0,1,0,1)), ((0,0),)*4)

print(id.reverse())
print(id)
assert id.reverse() == id
assert double.reverse() == RowColMap({0:0, 2:1, 4:2, 6:3}, {0:0, 2:1, 4:2}, False)

assert double.is_mappable_cell((1, 2))
assert not double.is_mappable_cell((1, 4))
assert not double.is_mappable_cell((4, 1))

assert double.is_mappable_gc(GriddedChord(Chord((0, 1, 1, 0)), ((0,0), (1,1), (1,1), (2,0))))
assert not double.is_mappable_gc(GriddedChord(Chord((0, 1, 1, 0)), ((0,0), (1, 1), (1, 1), (4,0))))
assert not double.is_mappable_gc(GriddedChord(Chord((1, 1)), ((0,0), (6,0))))

assert id.map_gc(crossed) == crossed