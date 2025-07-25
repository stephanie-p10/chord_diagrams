import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from strategies.requirement_insertion import RequirementInsertionStrategy, RequirementInsertionFactory
from chords import GriddedChord, Chord
from tiling import Tiling

crossed_chord = GriddedChord(Chord((0, 1, 0, 1)), ((0, 0),) *4)
def single_chord(pos_left, pos_right):
    return GriddedChord(Chord((0,0)), (pos_left, pos_right))

non_crossing = Tiling(obstructions=(crossed_chord,))
tiling1 = Tiling((single_chord((0,0),(2,0)), 
             single_chord((1,0), (2,0)), 
             single_chord((2,0), (2,0)), 
             GriddedChord(Chord((0,1)), ((2,0), (2,0))), 
             GriddedChord(Chord((1,0)), ((0,1), (1,1))), 
             GriddedChord(Chord((1,0)), ((2,0), (2,0))),), 
             ((single_chord((0,0), (1,0)),), 
              (single_chord((1,1), (2,1)),)))
tiling2 = Tiling((GriddedChord(Chord((0,1,0,1)), ((1,1),) *4),), 
                 ((single_chord((0,0), (2,0)),), (GriddedChord(Chord((0,1)), ((0,0), (1,1))),)))
tilings = [non_crossing, tiling1, tiling2]

insert_single_chord = RequirementInsertionStrategy((single_chord((0,0), (0,0)),))
decomp_nc = insert_single_chord.decomposition_function(non_crossing)
decomp_tiling1 = insert_single_chord.decomposition_function(tiling1)
decomp_tiling2 = insert_single_chord.decomposition_function(tiling2)
decomp_list = [decomp_nc, decomp_tiling1, decomp_tiling2]

def print_decomp(decomp):
    for child in decomp:
        print(child)

#for decomp in decomp_list:
#    print_decomp(decomp)
#    print()

obj_nc_contains = GriddedChord(Chord((0,0,1,1)), ((0,0),) * 4)
obj_nc_avoids = GriddedChord.empty_chord()

obj_t1_contains = GriddedChord(Chord((0,0,1,2,1,2)), ((0,0), (0,0), (0,0), (1,1), (1,0), (2,1)))
obj_t1_avoids = GriddedChord(Chord((0,1,0,1)), ((0,0), (1,1), (1,0), (2,1)))

obj_t2_contains = GriddedChord(Chord((0,1,0,2,2,1)), ((0,0), (0,0), (0,0), (1, 1), (1, 1), (2, 0)))
obj_t2_avoids = GriddedChord(Chord((0, 1, 1, 0)), ((0,0), (1, 1), (1, 1), (2, 0)))

def print_fwd(cls, obj):
    mapped_obj = insert_single_chord.forward_map(cls, obj)
    for obj in mapped_obj:
        print(obj)
    print()

obj_list = [obj_nc_avoids, obj_nc_contains, obj_t1_avoids, obj_t1_contains, obj_t2_avoids, obj_t2_contains]

for i in range(6):
    print_fwd(tilings[i//2], obj_list[i])

print_fwd(non_crossing, obj_nc_avoids)

assert insert_single_chord.backward_map(non_crossing, insert_single_chord.forward_map(non_crossing, obj_nc_contains)) == obj_nc_contains
assert insert_single_chord.backward_map(non_crossing, insert_single_chord.forward_map(non_crossing, obj_nc_avoids)) == obj_nc_avoids


req_ins_factory = RequirementInsertionFactory()

# sToDo: put in assert
print("Testing factory")
for strat in req_ins_factory(non_crossing):
    print("Strategy:")
    for part in strat.decomposition_function(non_crossing):
        print(part)
        print()
    print()

