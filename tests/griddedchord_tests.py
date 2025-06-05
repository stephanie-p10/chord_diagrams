import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from chords import *

# GriddedChord tests
gc1 = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 3), (2, 4), (2, 1), (3, 3), (3, 1), (4, 4)))
gc2 = GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)))
gc3 = GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (1, 0), (2, 1)))
gc4 = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), ((0, 0), (1, 1), (1, 1), (1, 0), (2, 3), (2, 4), (2, 1), (3, 3), (2, 1), (4, 4)))
empty_c = GriddedChord()

assert(GriddedChord.single_cell(Chord((0, 1, 0, 2, 2, 1)), (1, 1)) == gc2)
assert(GriddedChord.single_cell(Chord((0, 0, 1, 1, 2, 2)), (0, 0)) == 
                                GriddedChord(Chord((0, 0, 1, 1, 2, 2)), ((0,0), (0,0), (0,0), (0,0), (0,0), (0,0))))

assert(GriddedChord.empty_chord() == empty_c)
assert(GriddedChord.empty_chord() == GriddedChord(Chord(), ()))

assert(GriddedChord.single_chord(((0,0),)) == GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))))

assert(GriddedChord.single_chord(((0, 1), (1, 1))) == GriddedChord(Chord((0, 0)), ((0, 1), (1, 1))))

assert(gc2.occupies((1, 1)))
assert(not gc2.occupies((0, 0)))
assert(not gc3.occupies((2, 3)))
assert(gc1.occupies((2, 3)))

assert(set(gc3.occurrences_in(gc1)) == set([(0, 2)]))
assert(set(gc3.occurrences_in(gc4)) == set([(0, 2), (0, 1)]))
assert(set(gc3.occurrences_in(gc2)) == set())

assert(gc3.occurs_in(gc1))
assert(gc3.occurs_in(gc4))
assert(not gc3.occurs_in(gc2))
assert(not gc1.occurs_in(gc3))

assert gc1.contains(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (3, 1))), 
                    GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (3, 1), (1, 0))))
assert gc1.contains(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (3, 1))))
assert not gc1.contains(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (3, 1), (1, 0))))

assert not gc1.avoids(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (3, 1))), 
                    GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (3, 1), (1, 0))))
assert not gc1.avoids(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (3, 1)))) 
assert gc1.avoids(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (3, 1), (1, 0))))

assert not GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (1,0), (2,1))).contradictory()
assert GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (0,1), (2,1))).contradictory()
assert GriddedChord(Chord((0, 1, 0, 1)), ((0,1), (1,0), (1,1), (2,0))).contradictory()
assert GriddedChord(Chord((0, 1, 0, 1)), ((0, 1), (1, 0), (1, 1), (2, 0))).contradictory()

assert type(gc1.remove_cells([(0, 0)])) == GriddedChord
assert gc1.remove_cells([(0, 0)]) == GriddedChord(Chord((0, 1, 2, 3, 1, 2, 0, 3)), ((0, 1), (1, 1), (2, 3), (2, 4), (2, 1), (3, 3), (3, 1), (4, 4)))
assert gc2.remove_cells([(1, 1)]) == GriddedChord(Chord(), ())
assert gc1.remove_cells([(0, 0), (1, 1)]) == GriddedChord(Chord((0, 1, 2, 1, 0, 2)), ((0, 1), (2, 3), (2, 4), (3, 3), (3, 1), (4, 4)))

assert list(gc1.points_in_cell((0, 0))) == [0]
assert list(gc1.points_in_cell((1, 0))) == [3]
assert list (gc2.points_in_cell((1, 1))) == [0, 1, 2, 3, 4, 5]

assert set(GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).isolated_cells()) == set([(0, 0)])
assert set(GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (2, 2), (2, 2), (3, 1))).isolated_cells()) == set([(0, 0), (1, 1), (2, 2), (3, 1)])
assert set(GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (2, 1), (2, 2), (2, 1))).isolated_cells()) == set([(0, 0)])

assert GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([0])
assert not GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([1, 0])
assert not GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([1])
assert not GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (2, 2), (2, 2), (2, 1))).is_isolated([3])

assert GriddedChord(Chord((0, 1, 0, 1, 2, 2)), ((0,0), (0,0), (0,0), (0,0), (0, 0), (0, 0))).forced_point_index((0, 0), DIR_EAST) == 2
assert GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0,0), (0,0), (0,0), (0,0), (0, 0), (0, 0))).forced_point_index((0, 0), DIR_NORTH) == 2
assert GriddedChord(Chord((0, 1, 0, 1, 2, 2)), ((0,0), (0,0), (0,0), (0,0), (0, 0), (0, 0))).forced_point_index((0, 0), DIR_WEST) == 0
assert GriddedChord(Chord((0, 1, 0, 1, 2, 2)), ((0,0), (0,0), (0,0), (0,0), (0, 0), (0, 0))).forced_point_index((0, 0), DIR_SOUTH) == 0

assert set(gc1.get_points_row(1)) == set([(1, 1), (8, 1), (2, 2), (6, 2)])
assert set(gc1.get_points_row(2)) == set()
assert set(gc1.get_points_row(3)) == set([(4, 3), (7, 3)])

assert set(gc1.get_points_col(0)) == set([(0, 0), (1, 1)])
assert set(gc1.get_points_col(5)) == set()
assert set(gc1.get_points_col(1)) == set([(3, 0), (2, 2)])

assert set(gc1.get_points_left_col(2)) == set([(0, 0), (3, 0), (1, 1), (2, 2)])
assert set(gc1.get_points_left_col(0)) == set()

assert set(gc1.get_points_right_col(2)) == set([(8, 1), (7, 3), (9, 4)])
assert set(gc1.get_points_right_col(4)) == set()

assert set(gc1.get_points_below_row(0)) == set()
assert set(gc1.get_points_below_row(1)) == set([(0, 0), (3, 0)])
assert set(gc1.get_points_below_row(2)) == set([(0, 0), (3, 0), (1, 1), (8, 1), (2, 2), (6, 2)])

assert set(gc1.get_points_above_row(4)) == set()
assert set(gc1.get_points_above_row(3)) == set([(5, 4), (9, 4)])
assert set(gc1.get_points_above_row(2)) == set([(5, 4), (9, 4), (4, 3), (7, 3)])

assert gc1.get_subchord_below_row(0) == GriddedChord()
assert gc1.get_subchord_below_row(1) == GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))
assert gc1.get_subchord_below_row(2) == GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (3, 1)))

assert gc1.get_subgrid_at_chords([0, 3]) == GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (1, 0), (2, 3), (3, 3)))
assert gc1.get_subgrid_at_chords([2, 3]) == GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (2, 3), (2, 1), (3, 3)))
assert gc1.get_subgrid_at_chords([]) == GriddedChord()

assert GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (3, 1))).remove_chord_idx(0) == GriddedChord(Chord((0, 1, 1, 0)), ((0, 1), (1, 1), (2, 1), (3, 1)))
assert GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (3, 1))).remove_chord_idx(3) == GriddedChord(Chord((0, 1, 1, 0)), ((0, 1), (1, 1), (2, 1), (3, 1)))
assert GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (3, 1))).remove_chord_idx(4) == GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (3, 1)))

assert gc1.get_subchord_in_cells(((0, 0), (0, 1), (1, 0))) == GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (3, 1)))
assert gc1.get_subchord_in_cells(((4, 4), (3, 3))) == GriddedChord(Chord((0, 1, 0, 1)), ((2, 3), (2, 4), (3, 3), (4, 4)))
assert gc1.get_subchord_in_cells([(2, 0)]) == GriddedChord()

assert list(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (2, 1))).get_bounding_box(0, 1)) == [0, 2, 2, 3]
assert list(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 1), (1, 0), (1, 1))).get_bounding_box(0, 1)) == [0, 1, 1, 4]
assert list(gc1.get_bounding_box(1, 2)) == [2, 4, 4, 7]
assert list(gc1.get_bounding_box(0, 5)) == [0, 2, 10, 10]
assert list(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 1), (1, 0), (2, 1))).get_bounding_box(0, 2)) == [0, 1, 3, 4]

assert GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (2, 0))).insert_specific_chord(1, 1, 1, 2, 2) == GriddedChord(Chord((0, 1, 2, 2, 1, 0)), ((0, 0), (0, 1), (1, 1), (1, 1), (1, 1), (2, 0)))
assert GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (2, 0))).insert_specific_chord(1, 0, 1, 1, 2) == GriddedChord(Chord((0, 1, 2, 1, 2, 0)), ((0, 0), (0, 1), (0, 1), (1, 1), (1, 1), (2, 0)))
assert GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (2, 0))).insert_specific_chord(3, 1, 1, 3, 3) == GriddedChord(Chord((0, 1, 1, 2, 2, 0)), ((0, 0), (0, 1), (1, 1), (1, 3), (1, 3), (2, 0)))
assert GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (2, 0))).insert_specific_chord(3, 2, 2, 2, 2) == None

assert set(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (1, 0))).insert_chord(0, 0, 1)) == set([
        GriddedChord(Chord((0, 1, 2, 2, 0, 1)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 2, 1, 0)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 1, 2, 0)), ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 2, 1, 0)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 2, 0, 1)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0)))])

assert GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (1, 0))).insert_chord(0, 1, 1) == []
assert set(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (1, 0))).insert_chord(1, 0, 0)) == set([
    GriddedChord(Chord((0, 1, 1, 2, 2, 0)), ((0, 0), (0, 1), (0, 1), (0, 1), (1, 1), (1, 0))),
    GriddedChord(Chord((0, 1, 2, 1, 2, 0)), ((0, 0), (0, 1), (0, 1), (0, 1), (1, 1), (1, 0))),
    GriddedChord(Chord((0, 1, 2, 2, 1, 0)), ((0, 0), (0, 1), (0, 1), (0, 1), (1, 1), (1, 0)))
    ])
assert GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (1, 0))).insert_chord(1, 2, 2) == [
    GriddedChord(Chord((0, 1, 1, 0, 2, 2)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (2, 1)))]
assert set(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (1, 0))).insert_chord(0, 0, 1)) == set([
    GriddedChord(Chord((0, 1, 2, 2, 0, 1)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
    GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (1, 0))),
    GriddedChord(Chord((0, 1, 2, 2, 1, 0)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
    GriddedChord(Chord((0, 1, 2, 1, 2, 0)), ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (1, 0))),
    GriddedChord(Chord((0, 1, 2, 2, 1, 0)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
    GriddedChord(Chord((0, 1, 2, 2, 0, 1)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0)))])

assert GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (3, 1))).remove_chord(0) == GriddedChord(Chord((0, 1, 1, 0)), ((0, 1), (1, 1), (2, 1), (3, 1)))
assert GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (3, 1))).remove_chord(3) == GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (3, 1)))
assert GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 1), (1, 1), (1, 0), (2, 1), (3, 1))).remove_chord(2) == GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (3, 1)))

assert set(gc3.all_subchords()) == set([GriddedChord(), GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))), GriddedChord(Chord((0, 0)), ((1, 1), (2, 1)))])
assert set(GriddedChord(Chord((0, 1, 2, 2, 0, 1)), ((0, 0), (1, 0), (1, 1), (1, 1), (2, 0), (2, 0))).all_subchords(proper=False)) == set([
    GriddedChord(),
    GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),
    GriddedChord(Chord((0, 0)), ((1, 0), (2, 0))),
    GriddedChord(Chord((0, 0)), ((1, 1), (1, 1))),
    GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (2, 0), (2, 0))),
    GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 1), (1, 1), (2, 0))),
    GriddedChord(Chord((0, 1, 1, 0)), ((1, 0), (1, 1), (1, 1), (2, 0))),
    GriddedChord(Chord((0, 1, 2, 2, 0, 1)), ((0, 0), (1, 0), (1, 1), (1, 1), (2, 0), (2, 0)))
])

def add_one_x(cell: Cell):
    return (cell[0] + 1, cell[1])
assert gc3.apply_map(add_one_x) == GriddedChord(Chord((0, 1, 0, 1)), ((1,0), (2,1), (2, 0), (3, 1)))

assert not gc3.is_single_chord()
assert GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))).is_single_chord()

assert not gc1.is_single_cell()
assert gc2.is_single_cell()

assert not gc1.is_single_row()
assert gc2.is_single_row()
assert GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (2, 1))).is_single_row()

assert not gc1.is_single_col()
assert gc2.is_single_col()
assert GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 2), (1, 2), (1, 1))).is_single_col()

assert empty_c.is_empty()
assert not gc1.is_empty()

assert not gc2.is_interleaving()
#assert not GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 2))).is_interleaving()
assert GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 2), (1, 2), (2, 1))).is_interleaving()

assert gc1.compress() == [0, 1, 2, 0, 3, 4, 2, 3, 1, 4, 0, 0, 0, 1, 1, 1, 1, 0, 2, 3, 2, 4, 2, 1, 3, 3, 3, 1, 4, 4]
assert GriddedChord.decompress([0, 1, 2, 0, 3, 4, 2, 3, 1, 4, 0, 0, 0, 1, 1, 1, 1, 0, 2, 3, 2, 4, 2, 1, 3, 3, 3, 1, 4, 4]) == gc1
assert GriddedChord.decompress(gc2.compress()) == gc2

assert gc1.to_jsonable() == {'patt': (0, 1, 2, 0, 3, 4, 2, 3, 1, 4),
                             'pos': ((0, 0), (0, 1), (1, 1), (1, 0), (2, 3), (2, 4), (2, 1), (3, 3), (3, 1), (4, 4))}
assert GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (2, 1))).to_jsonable() == {
    'patt': (0, 1, 0, 2, 2, 1),
    'pos': ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (2, 1))
}

#assert GriddedChord.from_json("{'patt': (0, 1, 2, 0, 3, 4, 2, 3, 1, 4), 'pos': ((0, 0), (1, 0), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4))}") == gc1
assert GriddedChord.from_dict({'patt': (0, 1, 2, 0, 3, 4, 2, 3, 1, 4), 'pos': ((0, 0), (0, 1), (1, 1), (1, 0), (2, 3), (2, 4), (2, 1), (3, 3), (3, 1), (4, 4))}) == gc1
assert GriddedChord.from_dict(gc1.to_jsonable()) == gc1

assert gc1.patt == (0, 1, 2, 0, 3, 4, 2, 3, 1, 4)
assert gc1.pos == ((0, 0), (0, 1), (1, 1), (1, 0), (2, 3), (2, 4), (2, 1), (3, 3), (3, 1), (4, 4))

# sToDo: more tests on this
assert set(GriddedChord.all_grids(Chord((0, 0, 1, 1)), [(0, 0)])) == set([GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))])
assert set(GriddedChord.all_grids(Chord((0, 0, 1, 1)), [(0, 0), (0, 1)])) == set([
        GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 1), (0, 1))),
        GriddedChord(Chord((0, 0, 1, 1)), ((0, 1), (0, 1), (0, 1), (0, 1))),
        GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))])

print("all assertions passed")
