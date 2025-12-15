from tiling import Tiling, Chord, GriddedChord

empty_T = Tiling((GriddedChord(Chord((0,0)), ((0,0), (0,0))),), (), ())

print(empty_T.is_atom())
print(empty_T.minimum_size_of_object())