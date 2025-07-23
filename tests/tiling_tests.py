import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

import time

from tiling import Tiling
from chords import *

gc_empty = GriddedChord()

# For the following: gc_single_ab_cd is single chord that starts in (a, b), ends in (c, d)
gc_single_00_00 = GriddedChord(Chord((0, 0)), ((0, 0), (0, 0)))
gc_single_00_10 = GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))
gc_single_10_10 = GriddedChord(Chord((0, 0)), ((1, 0), (1, 0)))
gc_single_01_11 = GriddedChord(Chord((0, 0)), ((0, 1), (1, 1)))
gc_single_11_11 = GriddedChord(Chord((0, 0)), ((1, 1), (1, 1)))

# For the following: chord 0: (0, 0), (1, 0); chord 1: (0, 1), (1, 1) 
gc_crossed = GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (1, 1)))
gc_nested = GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 1), (1, 1), (1, 0)))
gc_3_crossed = gc_crossed.insert_specific_chord(2, 0, 1, 2, 4)
gc_3_nested = gc_nested.insert_specific_chord(2, 0, 1, 2, 2)

gc_crossed_h1 = GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0)))
gc_3_crossed_h1 = gc_crossed_h1.insert_specific_chord(0, 0, 1, 2, 4)
gc_nested_h1 = GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0)))
gc_3_nested_h1 = gc_nested_h1.insert_specific_chord(0, 0, 1, 2, 2)

gc_disjoint = GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (1, 1), (1, 1)))

# For the following single cell chords, the cell is (0,0)
gc_sc_crossed = GriddedChord.single_cell(Chord((0, 1, 0, 1)), (0, 0))
gc_sc_nested = GriddedChord.single_cell(Chord((0, 1, 1, 0)), (0, 0))
gc_sc_disjoint = GriddedChord.single_cell(Chord((0, 0, 1, 1)), (0, 0))

# Tilings:
t_no_restrictions = Tiling()
t_av_single = Tiling((gc_single_00_00,), (), ())
t_cn_single = Tiling((), ((gc_single_00_00,),))
t_lk_single = Tiling( (), (), (((0, 0),),))
t_lk_2x2 = Tiling((), (), (((0, 0), (0, 1), (1, 0), (1, 1)),))
t_ex = Tiling((gc_single_00_00, gc_single_11_11), ((gc_single_01_11, gc_disjoint), (gc_nested,)), (((0, 0), (0, 1), (1, 0)), ((1, 0), (1, 1))))

# contains tests
assert t_no_restrictions.contains(gc_single_00_00)
assert t_cn_single.contains(gc_single_00_00)
assert t_cn_single.contains(gc_sc_crossed)
assert not Tiling( (gc_nested,), ((gc_single_00_00,),)).contains(gc_crossed)
assert Tiling( (gc_nested,), ((gc_single_00_10,),)).contains(gc_crossed)
assert Tiling((gc_3_nested,), ((gc_single_00_10,),)).contains(gc_nested)
assert not Tiling((gc_nested,), ((gc_single_00_10,),)).contains(gc_3_nested)
assert not t_lk_single.contains(gc_single_00_10)
assert t_lk_single.contains(gc_single_00_00)
assert not t_lk_single.contains(gc_sc_disjoint)
assert t_lk_2x2.contains(gc_crossed)
assert not t_lk_2x2.contains(gc_3_nested)

# is_empty tests
#assert Tiling((1, 1), (gc_empty,), (), ()).is_empty()
#assert Tiling((1, 1), (gc_single_00_00,), (), ()).is_empty()
#assert Tiling((1, 1), (gc_sc_crossed,), ((gc_sc_crossed, gc_sc_disjoint, gc_sc_nested),), (((0, 0),),)).is_empty()
#assert Tiling((2, 1), (gc_single_00_00,), (), (((0, 0),),)).is_empty()
#assert Tiling((3, 2), (gc_crossed,), ((gc_3_crossed,),), ()).is_empty() # uh oh... this test takes a long time...
#assert Tiling((2, 1), (gc_crossed_h1,), ((gc_3_crossed_h1,),), ()).is_empty() 
#assert Tiling((2, 1), (gc_nested_h1,), ((gc_3_nested_h1,),), ()).is_empty()
#assert Tiling((2, 1), (gc_single_00_10, gc_single_00_00, gc_single_10_10), (), (((0, 0), (1, 0)),)).is_empty()
#assert Tiling((2, 1), (gc_single_00_00,), (), (((0, 0),),)).is_empty()

#assert not t_no_restrictions.is_empty()
#assert not Tiling((2, 2), (), (), ()).is_empty()
#assert not Tiling((2, 1), (gc_single_00_10,), (), (((0, 0), (1, 0)),)).is_empty()
#assert not Tiling((2, 1), (gc_3_nested_h1,), ((gc_nested_h1,),), ()).is_empty()
#assert Tiling((1, 1), (), ((gc_crossed, gc_nested), (gc_nested,)), ()).is_empty()
#start1 = time.time()
#assert not Tiling((2, 2), (gc_3_nested_h1,), ((gc_nested_h1,), (gc_crossed, gc_nested)), ()).is_empty_1()
#end1 = time.time()
#start2 = time.time()
#assert not Tiling((2, 2), (gc_3_nested_h1,), ((gc_nested_h1,), (gc_crossed, gc_nested)), ()).is_empty_2()
#end2 = time.time()

#assert Tiling.from_dict(t_ex.to_jsonable()) == t_ex
#print("is_empty_1: ", end1-start1)
#print("is_empty_2: ", end2-start2)

t = Tiling((GriddedChord(Chord((0,1)), ((0,1), (0,1))), GriddedChord(Chord((0,0)), ((0,1), (0,1))), GriddedChord(Chord((1,0)), ((0,1), (0,1))), GriddedChord(Chord((0,1)), ((0,0), (0,0))), GriddedChord(Chord((0,0)), ((0,0), (0,0))), GriddedChord(Chord((1,0)), ((0,0), (0,0)))),
           ((GriddedChord(Chord((0,)), ((0,0),)),),(GriddedChord(Chord((0,)), ((0,1),)),)))

small_tiling = Tiling((),
                      ((GriddedChord(Chord((0,1)), ((0,0), (1,1))),),))

large_list = small_tiling.all_chords_on_tiling(3, True)



#for chord in large_list:
    #print(chord)

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

#print()
#print(non_crossing.contains(GriddedChord(Chord((1, 0)), ((3, 1), (3, 1)))))

for chord in non_crossing.all_chords_on_tiling(4, False):
    print(chord)

print()
print(non_crossing.contains(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 1), (1, 1), (2, 0)))))
#print(t.point_cells)

print("asserts passed")




