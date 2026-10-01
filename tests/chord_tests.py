import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))
from src.common.chords import Chord

c1 = Chord((0, 1, 1, 0))
c2 = Chord((0, 1, 0, 1))
c3 = Chord((0, 1, 2, 2, 0, 1))

# Chord tests:

# Testing whether the length of the chord is correct
assert(len(c1) == 2)
assert(len(c3) == 3)
assert(len(Chord())== 0)

# Testing whether the pattern of the chord is correct
assert(c1.get_pattern() == ((0, 1, 1, 0)))
assert(c2.get_pattern() == (0, 1, 0, 1))

# Testing whether the concatenation of the chord is correct
assert(c1.concat(c2) == Chord((0, 1, 1, 0, 2, 3, 2, 3)))
assert(c1.concat(Chord()) == c1)
assert(Chord((0, 0, 1, 1)).concat(Chord((0, 1, 1, 0))) == Chord((0, 0, 1, 1, 2, 3, 3, 2)))

# Testing whether the removal of the chord is correct
#   (the input is the label of the chord to be removed)
assert(Chord((0, 1, 0, 1, 2, 2)).remove_chord()==Chord((0, 1, 0, 1)))
assert(Chord((0, 1, 1, 0, 2, 2)).remove_chord(2)==Chord((0, 1, 1, 0)))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2)).remove_chord(3)==Chord((0, 1, 1, 0, 2, 2)))

# Testing whether the removal of the index of the chord is correct 
#   (index is the index of the chord in the list of chords)
assert(Chord((0, 1, 0, 1, 2, 2)).remove_index()==Chord((0, 1, 0, 1)))
assert(Chord((0, 1, 1, 0, 2, 2)).remove_index(2)==Chord((0, 0, 1, 1)))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2)).remove_index(3)==Chord((0, 0, 1, 2, 2, 1)))

# Testing whether the conversion to standard chord is correct
assert(Chord.to_standard(("abba")) == c1)
assert(Chord.to_standard((2, -1, 2, -1)) == c2)
assert(Chord.to_standard((0, 4, 7, 7, 0, 4)) == c3)

# Testing whether the contains method is correct
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2))._contains(Chord((0, 1, 1, 0))))
assert(not Chord((0, 1, 1, 0, 2, 3, 3, 2))._contains(Chord((0, 1, 0, 1))))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2))._contains(Chord((0, 0, 1, 1))))

assert(not Chord((0,0,1,1)).contains(Chord((0,0,1,1)), Chord((0,1,0,1))))
assert(Chord((0, 1, 1, 0, 2, 3, 2, 3)).contains(Chord((0, 0, 1, 1)), Chord((0, 1, 1, 0))))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2)).contains(Chord((0, 1, 1, 0))))

# Testing whether the occurrences in method is correct
assert(set(Chord((0, 1, 0, 1)).occurrences_in(Chord((0, 1, 2, 0, 1, 2)))) == set([(0, 1), (1, 2), (0, 2)]))
assert(set(Chord((0, 0)).occurrences_in(Chord((0, 1, 1, 0)))) == set([(0,), (1,)]))
assert(set(Chord((0, 1, 1, 0)).occurrences_in(Chord((0, 1, 1, 0, 2, 3, 3, 2)))) == set([(0, 1), (2, 3)]))

# Testing whether the avoids method is correct
assert(not Chord((0, 0, 1, 1)).avoids(Chord((0, 0, 1, 1)), Chord((0, 1, 0, 1))))
assert(Chord((0, 1, 1, 0, 2, 3, 3, 2)).avoids(Chord((0, 1, 0, 1)), Chord((0, 1, 2, 2, 1, 0))))
assert(not Chord((0, 1, 1, 2, 0, 3, 3, 2)).avoids(Chord((0, 1, 1, 0))))

# Testing whether the insert method is correct
assert(Chord((0, 1, 0, 1)).insert() == Chord((0, 1, 0, 1, 2, 2)))
assert(Chord((0, 1, 0, 1)).insert(0, 0) == Chord((0, 0, 1, 2, 1, 2)))
assert(Chord((0, 1, 2, 2, 0, 1)).insert(0, 4) == Chord((0, 1, 2, 3, 3, 0, 1, 2)))
assert(Chord((0, 1, 2, 2, 0, 1)).insert(4, 4) == Chord((0, 1, 2, 2, 3, 3, 0, 1)))
assert(Chord((0, 1, 2, 2, 0, 1)).insert(4, 5) == Chord((0, 1, 2, 2, 3, 0, 3, 1)))

# Testing whether the get_chords method is correct
assert(Chord((0, 0, 1, 1)).get_chords([0]) == Chord((0, 0)))
assert(Chord((0, 1, 0, 1)).get_chords([0, 1]) == Chord((0, 1, 0, 1)))
assert(Chord((0, 1, 1, 0, 2, 2)).get_chords() == Chord((0, 1, 1, 0, 2, 2)))
assert(Chord((0, 1, 1, 0, 2, 2)).get_chords([1, 2]) == Chord((0, 0, 1, 1)))

# Testing whether the chord_dict method is correct
assert(Chord((0, 1, 1, 0, 2, 2)).chord_dict == {0: (0, 3), 1: (1, 2), 2: (4, 5)})

# Testing whether the connected method is correct
assert Chord((0, 1, 2, 3, 1, 3, 0, 4, 5, 2, 5, 4)).connected()
assert not Chord((0, 1, 2, 1, 2, 0, 3, 4, 4, 3)).connected()
assert not Chord.to_standard((0, 1, 2, 1, 2, 0, 3, 4, 4, 3)).connected()
assert Chord((0, 0)).connected()
assert not Chord().connected()

print("tests passed")