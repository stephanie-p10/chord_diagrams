from chords import *

c1 = Chord((0, 1, 1, 0))
c2 = Chord((0, 1, 0, 1))
c3 = Chord((0, 1, 2, 2, 0, 1))

# Chord tests:
assert(len(c1) == 2)
assert(len(c3) == 3)
assert(len(Chord())== 0)

assert(c1.get_pattern() == ((0, 1, 1, 0)))
assert(c2.get_pattern() == (0, 1, 0, 1))

assert(c1.concat(c2) == Chord((0, 1, 1, 0, 2, 3, 2, 3)))
assert(c1.concat(Chord()) == c1)
assert(Chord((0, 0, 1, 1)).concat(Chord((0, 1, 1, 0))) == Chord((0, 0, 1, 1, 2, 3, 3, 2)))

assert(Chord((0, 1, 0, 1, 2, 2)).remove_chord()==Chord((0, 1, 0, 1)))
assert(Chord((0, 1, 1, 0, 2, 2)).remove_chord(2)==Chord((0, 1, 1, 0)))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2)).remove_chord(3)==Chord((0, 1, 1, 0, 2, 2)))

assert(Chord((0, 1, 0, 1, 2, 2)).remove_index()==Chord((0, 1, 0, 1)))
assert(Chord((0, 1, 1, 0, 2, 2)).remove_index(2)==Chord((0, 0, 1, 1)))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2)).remove_index(3)==Chord((0, 0, 1, 2, 2, 1)))

assert(Chord.to_standard(("abba")) == c1)
assert(Chord.to_standard((2, -1, 2, -1)) == c2)
assert(Chord.to_standard((0, 4, 7, 7, 0, 4)) == c3)

assert(Chord((0, 1, 1, 0, 2, 3, 3, 2))._contains(Chord((0, 1, 1, 0))))
assert(not Chord((0, 1, 1, 0, 2, 3, 3, 2))._contains(Chord((0, 1, 0, 1))))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2))._contains(Chord((0, 0, 1, 1))))

assert(not Chord((0,0,1,1)).contains(Chord((0,0,1,1)), Chord((0,1,0,1))))
assert(Chord((0, 1, 1, 0, 2, 3, 2, 3)).contains(Chord((0, 0, 1, 1)), Chord((0, 1, 1, 0))))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2)).contains(Chord((0, 1, 1, 0))))

assert(set(Chord((0, 1, 0, 1)).occurrences_in(Chord((0, 1, 2, 0, 1, 2)))) == set([(0, 1), (1, 2), (0, 2)]))
assert(set(Chord((0, 0)).occurrences_in(Chord((0, 1, 1, 0)))) == set([(0,), (1,)]))
assert(set(Chord((0, 1, 1, 0)).occurrences_in(Chord((0, 1, 1, 0, 2, 3, 3, 2)))) == set([(0, 1), (2, 3)]))

assert(not Chord((0, 0, 1, 1)).avoids(Chord((0, 0, 1, 1)), Chord((0, 1, 0, 1))))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2)).avoids(Chord((0, 1, 0, 1)), Chord((0, 1, 2, 2, 1, 0))))
assert(not Chord((0, 1, 1, 2, 0, 3, 3, 2)).avoids(Chord((0, 1, 1, 0))))

assert(Chord((0, 1, 0, 1)).insert() == Chord((0, 1, 0, 1, 2, 2)))
assert(Chord((0, 1, 0, 1)).insert(0, 0) == Chord((0, 0, 1, 2, 1, 2)))
assert(Chord((0, 1, 2, 2, 0, 1)).insert(0, 4) == Chord((0, 1, 2, 3, 3, 0, 1, 2)))
assert(Chord((0, 1, 2, 2, 0, 1)).insert(4, 4) == Chord((0, 1, 2, 2, 3, 3, 0, 1)))
assert(Chord((0, 1, 2, 2, 0, 1)).insert(4, 5) == Chord((0, 1, 2, 2, 3, 0, 3, 1)))

assert(Chord((0, 0, 1, 1)).get_chords([0]) == Chord((0, 0)))
assert(Chord((0, 1, 0, 1)).get_chords([0, 1]) == Chord((0, 1, 0, 1)))
assert(Chord((0, 1, 1, 0, 2, 2)).get_chords() == Chord((0, 1, 1, 0, 2, 2)))
assert(Chord((0, 1, 1, 0, 2, 2)).get_chords([1, 2]) == Chord((0, 0, 1, 1)))

assert(Chord((0, 1, 1, 0, 2, 2)).chord_dict == {0: (0, 3), 1: (1, 2), 2: (4, 5)})

# GriddedChord tests
gc1 = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), ((0, 0), (1, 0), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4)))
gc2 = GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)))
gc3 = GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (0,1), (1,2)))
gc4 = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), ((0, 0), (1, 1), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 2), (4, 4)))
empty_c = GriddedChord()

assert(GriddedChord.single_cell(Chord((0, 1, 0, 2, 2, 1)), (1, 1)) == gc2)
assert(GriddedChord.single_cell(Chord((0, 0, 1, 1, 2, 2)), (0,0)) == 
                                GriddedChord(Chord((0, 0, 1, 1, 2, 2)), ((0,0), (0,0), (0,0), (0,0), (0,0), (0,0))))

assert(GriddedChord.empty_chord() == empty_c)
assert(GriddedChord.empty_chord() == GriddedChord(Chord(), ()))

assert(GriddedChord.single_chord(((0,0),)) == GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))))
assert(GriddedChord.single_chord(((1, 0), (1, 1))) == GriddedChord(Chord((0, 0)), ((1, 0), (1, 1))))

assert(gc2.occupies((1, 1)))
assert(not gc2.occupies((0, 0)))
assert(not gc3.occupies((3, 2)))
assert(gc1.occupies((3, 2)))

assert(set(gc3.occurrences_in(gc1)) == set([(0, 2)]))
assert(set(gc3.occurrences_in(gc4)) == set([(0, 2), (0, 1)]))
assert(set(gc3.occurrences_in(gc2)) == set())

assert(gc3.occurs_in(gc1))
assert(gc3.occurs_in(gc4))
assert(not gc3.occurs_in(gc2))
assert(not gc1.occurs_in(gc3))

assert(gc1.contains(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3))), 
                    GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 3), (0, 1)))))
assert(gc1.contains(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3)))))
assert not gc1.contains(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 3), (0, 1))))

assert not gc1.avoids(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3))), 
                    GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 3), (0, 1))))
assert not gc1.avoids(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3)))) 
assert gc1.avoids(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 3), (0, 1))))

assert not GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (0,1), (1,2))).contradictory()
assert GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (1,0), (1,2))).contradictory()
assert GriddedChord(Chord((0, 1, 0, 1)), ((1,0), (0,1), (1,1), (0,2))).contradictory()
assert GriddedChord(Chord((0, 1, 0, 1)), ((1, 0), (0, 1), (1, 1), (1, 2))).contradictory()

assert type(gc1.remove_cells([(0, 0)])) == GriddedChord
assert gc1.remove_cells([(0, 0)]) == GriddedChord(Chord((0, 1, 2, 3, 1, 2, 0, 3)), ((1, 0), (1, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4)))
assert gc2.remove_cells([(1, 1)]) == GriddedChord(Chord(), ())
assert gc1.remove_cells([(0, 0), (1, 1)]) == GriddedChord(Chord((0, 1, 2, 1, 0, 2)), ((1, 0), (3, 2), (4, 2), (3, 3), (1, 3), (4, 4)))

assert list(gc1.points_in_cell((0, 0))) == [0]
assert list(gc1.points_in_cell((0, 1))) == [3]
assert list (gc2.points_in_cell((1, 1))) == [0, 1, 2, 3, 4, 5]

assert set(GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).isolated_cells()) == set([(0, 0)])
assert set(GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (2, 2), (2, 2), (1, 3))).isolated_cells()) == set([(0, 0), (1, 1), (2, 2), (1, 3)])
assert set(GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 2), (2, 2), (1, 2))).isolated_cells()) == set([(0, 0)])

assert GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([0])
assert not GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([0, 1])
assert not GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([1])
assert not GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (2, 2), (2, 2), (1, 2))).is_isolated([3])




print("all assertions passed")




