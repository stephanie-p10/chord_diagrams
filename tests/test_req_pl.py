import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from itertools import chain

import pytest

from misc import DIR_EAST, DIR_NORTH, DIR_SOUTH, DIR_WEST
from chords import GriddedChord, Chord
from tiling import Tiling
from algorithms.requirement_placement import RequirementPlacement

non_crossing = Tiling(
        obstructions=(
            GriddedChord(Chord((0,1,0,1)), ((0, 0),) * 4),
        ),
        requirements=(
            [GriddedChord(Chord((0,0)), ((0,0), (0,0)))],
        )
    )


tiling1 = Tiling(
        obstructions=[
            GriddedChord(Chord((0, 0)), ((1, 1), (1, 1))),
            GriddedChord(Chord((0, 1)), ((1, 1), (1, 1))),
            GriddedChord(Chord((1, 0)), ((1, 1), (1, 1))),
            GriddedChord(Chord((0, 0)), ((2, 1), (2, 1))),
            GriddedChord(Chord((0, 1)), ((2, 1), (2, 1))),
            GriddedChord(Chord((1, 0)), ((2, 1), (2, 1))),
            GriddedChord(Chord((0,)), ((1, 0),)),
            GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))),
            GriddedChord(Chord((0,)), ((0, 1),)),
        ],
        requirements=[
            [GriddedChord(Chord((0,)), ((1, 1),))],
        ],
    )

'''tiling2 = Tiling(
        obstructions=[
            GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), ) * 4),
            GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
            GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))
        ],
        requirements=[
            [GriddedChord(Chord((0, 0)), ((1, 1), (2, 1)))],
            [GriddedChord(Chord((0, 1)), ((0, 0), (1, 1)))],
        ],
    )'''
tiling3 = Tiling(
        obstructions=[GriddedChord(Chord((0, 1, 1, 0)), ((1, 1),) * 4)],
        requirements=[[GriddedChord(Chord((0,)), ((0, 0),))]],
    )

'''tiling4 = Tiling(obstructions=[GriddedChord(Chord((0,1)), ((0,1), (0,1))), 
                               GriddedChord(Chord((0,0)), ((0,1), (0,1))), 
                               GriddedChord(Chord((1,0)), ((0,1), (0,1))), 
                               GriddedChord(Chord((0,1)), ((0,0), (0,0))),
                               GriddedChord(Chord((0,0)), ((0,0), (0,0))), 
                               GriddedChord(Chord((1,0)), ((0,0), (0,0)))],
                requirements=[[GriddedChord(Chord((0,)), ((0,0),))],
                              [GriddedChord(Chord((0,)), ((0,1),))],
                              [GriddedChord(Chord((0,)), ((1,0),))],
                              [GriddedChord(Chord((0,)), ((1,1),))]])
'''
tiling5 = Tiling(obstructions=[GriddedChord(Chord((0, 1, 1, 0)), ((1, 1),) * 4),  
                               GriddedChord(Chord((0,1)), ((0,0), (0,0))),
                               GriddedChord(Chord((0,0)), ((0,0),(0,0))), 
                               GriddedChord(Chord((1,0)), ((0,0), (0,0)))],
                requirements=[[GriddedChord(Chord((0,)), ((0,0),))]])
print("done")

placement1 = RequirementPlacement(tiling1)
placement1ownrow = RequirementPlacement(tiling1, own_row=True, own_col=False)
placement1owncol = RequirementPlacement(tiling1, own_row=False, own_col=True)

#placement2 = RequirementPlacement(tiling2)
#placement2ownrow = RequirementPlacement(tiling2, own_row=True, own_col=False)
#placement2owncol = RequirementPlacement(tiling2, own_row=False, own_col=True)

placement3 = RequirementPlacement(tiling3)

#placement4 = RequirementPlacement(tiling4)

placement5 = RequirementPlacement(tiling5)

placement_nc = RequirementPlacement(non_crossing)
#placement_nc_ownroq = RequirementPlacement(tiling2, own_row=True, own_col=False)
#placement_nc_owncol = RequirementPlacement(tiling2, own_row=False, own_col=True)

gc1 = GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), (0, 0), (1, 0), (1, 1), (2, 0), (2,1)))

#assert placement1._tiling_point_col_cells() == frozenset([(1, 1)])
#assert placement2._tiling_point_col_cells() == frozenset([])
#assert placement1._tiling_point_row_cells() == frozenset([(1, 1), (2, 1)])
#assert placement5._tiling_point_row_cells() == frozenset([(0, 0)])

#assert placement1.already_placed(GriddedChord(Chord((0,0)), ((1,1), (2,1))), 0)
#assert not placement1.already_placed(GriddedChord(Chord((0,0)), ((1,1), (2,1))), 1)

#assert placement1._gridded_chord_translation(gc1, (2,1)) == GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), (0, 2), (3, 2), (3, 3), (4, 0), (4,3)))
#assert placement2._gridded_chord_translation(gc1, (3,2)) == GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), (0, 0), (1, 0), (3, 3), (4, 0), (4,3)))
#assert placement1ownrow._gridded_chord_translation(gc1, (2, 1)) == GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), (0, 2), (1, 2), (1, 3), (2, 0), (2,3)))
#assert placement1owncol._gridded_chord_translation(gc1, (2, 1)) == GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), (0, 0), (3, 0), (3, 1), (4, 0), (4,1)))

assert set(placement1.get_multiplexes_of_chord(GriddedChord(Chord((0,0)), ((0,0), (0,0))), (0,0))) == set([GriddedChord(Chord((0,0)), ((0,0), (0,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((2,2), (2,2))),
                                                                                                           GriddedChord(Chord((0,0)), ((2,0), (2,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((1,1), (2,1))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,1), (1,1))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,2), (0,2))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,0), (2,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,2), (2,2)))])
assert set(placement1ownrow.get_multiplexes_of_chord(GriddedChord(Chord((0,0)), ((0,0), (0,0))), (0,0))) == set([GriddedChord(Chord((0,0)), ((0,0), (0,0))),
                                                                                                                 GriddedChord(Chord((0,0)), ((0,2), (0,2))),
                                                                                                                 GriddedChord(Chord((0,0)), ((0,1), (0,1)))])
print("finished")
assert set(placement1owncol.get_multiplexes_of_chord(GriddedChord(Chord((0,0)), ((0,0), (0,0))), (0,0))) == set([GriddedChord(Chord((0,0)), ((0,0), (0,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((2,0), (2,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,0), (2,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((1,0), (2,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,0), (1,0)))])

assert set(placement1.get_multiplexes_of_chords([GriddedChord(Chord((0,0)), ((0,0), (0,0))), 
                                                GriddedChord(Chord((0,1)), ((0,0), (0,1)))], (0,0))) == set([GriddedChord(Chord((0,0)), ((0,0), (0,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((2,2), (2,2))),
                                                                                                           GriddedChord(Chord((0,0)), ((2,0), (2,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((1,1), (2,1))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,1), (1,1))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,2), (0,2))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,0), (2,0))),
                                                                                                           GriddedChord(Chord((0,0)), ((0,2), (2,2))),
                                                                                                           GriddedChord(Chord((0,1)), ((2,2), (2,3))),
                                                                                                           GriddedChord(Chord((0,1)), ((2,0), (2,3))),
                                                                                                           GriddedChord(Chord((0,1)), ((1,1), (2,3))),
                                                                                                           GriddedChord(Chord((0,1)), ((0,2), (2,3))),
                                                                                                           GriddedChord(Chord((0,1)), ((0,0), (2,3))),
                                                                                                           GriddedChord(Chord((0,1)), ((0,0), (0,3))),
                                                                                                           GriddedChord(Chord((0,1)), ((0,2), (0,3)))])

assert set(placement_nc.stretched_obs((0, 0),)) == set([GriddedChord(Chord((0, 1, 0, 1)), ((2, 2),) * 4),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 2),) * 4),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((2, 0),) * 4),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0),) * 4), 
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 2), (2, 2), (2, 2), (2, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 2), (0, 2), (2, 2), (2, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 2), (0, 2), (0, 2), (2, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (2, 0), (2, 0), (2, 0))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (2, 0))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (2, 2), (2, 1), (2, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 1), (2, 0), (2, 1))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 1), (0, 2), (1, 1), (2, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (0, 0), (1, 1))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (2, 2), (2, 0), (2, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 2), (2, 0), (2, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 2), (0, 0), (2, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 2), (0, 0), (0, 2))),
                                                      GriddedChord(Chord((0, 1, 0, 1)), ((2, 0), (2, 2), (2, 0), (2, 2)))])
'''assert set(placement1.stretched_obs((2, 0),)) == set([GriddedChord(Chord((0,)), ((0, 3),)),
                                                   GriddedChord(Chord((0,)), ((1, 0),)),
                                                   GriddedChord(Chord((0,)), ((1, 2),)),
                                                   GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))),
                                                   GriddedChord(Chord((0, 0)), ((1, 3), (1, 3))), 
                                                   GriddedChord(Chord((1, 0)), ((1, 3), (1, 3))), 
                                                   GriddedChord(Chord((0, 1)), ((1, 3), (1, 3))),
                                                   GriddedChord(Chord((1, 0)), ((0, 2), (0, 0))),
                                                   GriddedChord(Chord((1, 0)), ((0, 2), (0, 2))),

                                                   GriddedChord(Chord((0, 0)), ((2, 3), (2, 3))), 
                                                   GriddedChord(Chord((1, 0)), ((2, 3), (2, 3))), 
                                                   GriddedChord(Chord((0, 1)), ((2, 3), (2, 3))), 
                                                   
                                                   GriddedChord(Chord((0, 0)), ((2, 3), (4, 3))), 
                                                   GriddedChord(Chord((1, 0)), ((2, 3), (4, 3))), 
                                                   GriddedChord(Chord((0, 1)), ((2, 3), (4, 3))), 

                                                   GriddedChord(Chord((0, 0)), ((4, 3), (4, 3))), 
                                                   GriddedChord(Chord((1, 0)), ((4, 3), (4, 3))), 
                                                   GriddedChord(Chord((0, 1)), ((4, 3), (4, 3)))])'''

assert set(placement_nc.added_obs((0, 0), (0, 0), True)) == set([GriddedChord(Chord((0,)), ((0, 1),)),
                                                                 GriddedChord(Chord((0,)), ((1, 0),)),
                                                                 GriddedChord(Chord((0,)), ((1, 2),)),
                                                                                
                                                                 GriddedChord(Chord((0, 0)), ((1, 1), (1, 1))), 
                                                                 GriddedChord(Chord((1, 0)), ((1, 1), (1, 1))), 
                                                                 GriddedChord(Chord((0, 1)), ((1, 1), (1, 1))), 

                                                                 GriddedChord(Chord((0, 0)), ((2, 1), (2, 1))), 
                                                                 GriddedChord(Chord((1, 0)), ((2, 1), (2, 1))), 
                                                                 GriddedChord(Chord((0, 1)), ((2, 1), (2, 1)))])

assert set(placement_nc.added_obs((0, 0), (0, 0))) == set([GriddedChord(Chord((0,)), ((2, 1),)),
                                                           GriddedChord(Chord((0,)), ((1, 0),)),
                                                           GriddedChord(Chord((0,)), ((1, 2),)),
                                                                                
                                                           GriddedChord(Chord((0, 0)), ((1, 1), (1, 1))), 
                                                           GriddedChord(Chord((1, 0)), ((1, 1), (1, 1))), 
                                                           GriddedChord(Chord((0, 1)), ((1, 1), (1, 1))), 

                                                           GriddedChord(Chord((0, 0)), ((0, 1), (0, 1))), 
                                                           GriddedChord(Chord((1, 0)), ((0, 1), (0, 1))), 
                                                           GriddedChord(Chord((0, 1)), ((0, 1), (0, 1)))])
'''assert set(placement2.added_obs((1, 1), (2, 1), True)) == set([GriddedChord(Chord((0,)), ((0, 2),)),
                                                               GriddedChord(Chord((0,)), ((1, 2),)),
                                                               GriddedChord(Chord((0,)), ((3, 2),)),
                                                               GriddedChord(Chord((0,)), ((2, 1),)),
                                                               GriddedChord(Chord((0,)), ((2, 0),)),
                                                               GriddedChord(Chord((0,)), ((2, 3),)),
                                                                                
                                                               GriddedChord(Chord((0, 0)), ((2, 2), (2, 2))), 
                                                               GriddedChord(Chord((1, 0)), ((2, 2), (2, 2))), 
                                                               GriddedChord(Chord((0, 1)), ((2, 2), (2, 2))), 

                                                               GriddedChord(Chord((0, 0)), ((4, 2), (4, 2))), 
                                                               GriddedChord(Chord((1, 0)), ((4, 2), (4, 2))), 
                                                               GriddedChord(Chord((0, 1)), ((4, 2), (4, 2)))])'''
print("north: ", DIR_NORTH, ", south: ", DIR_SOUTH, ", east: ", DIR_EAST, ", west: ", DIR_WEST)
#print(placement2.point_multiplex_obs(gc1, True))

assert set(placement1.point_multiplex_obs(GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), )*6), 3)) == set(
    [GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((1, 1), (2, 2), (2, 2), (2, 2), (2, 1), (2, 2))), 
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (1, 1), (2, 1), (2, 2), (2, 0), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 1), (1, 1), (2, 2), (2, 0), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 1), (0, 2), (0, 2), (0, 2), (1, 1), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 0), (0, 0), (0, 1), (0, 0), (1, 1)))])

assert set(placement1.point_multiplex_obs(GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), )*6), 5)) == set(
    [GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((1, 1), (2, 2), (2, 2), (2, 2), (2, 1), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (1, 1), (2, 1), (2, 2), (2, 0), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 1), (1, 1), (2, 2), (2, 0), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 0), (0, 0), (1, 1), (2, 0), (2, 1))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 1), (0, 2), (0, 2), (0, 2), (1, 1), (2, 2)))])

assert set(placement1.point_multiplex_obs(GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), )*6), 4)) == set(
    [GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((1, 1), (2, 2), (2, 2), (2, 2), (2, 1), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (1, 1), (2, 1), (2, 2), (2, 0), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 1), (1, 1), (2, 2), (2, 0), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 0), (0, 0), (1, 1), (2, 0), (2, 1))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 0), (0, 0), (0, 1), (0, 0), (1, 1)))])

assert set(placement1.point_multiplex_obs(GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), )*6), 0)) == set(
    [GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 1), (0, 2), (0, 2), (0, 2), (1, 1), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (1, 1), (2, 1), (2, 2), (2, 0), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 1), (1, 1), (2, 2), (2, 0), (2, 2))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 0), (0, 0), (1, 1), (2, 0), (2, 1))),
        GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (0, 0), (0, 0), (0, 1), (0, 0), (1, 1)))])

'''assert placement2.point_multiplex_obs(gc1, 5) == []
assert placement2.point_multiplex_obs(gc1, 3) == []
assert placement2.point_multiplex_obs(gc1, 0) == [GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (1, 1), (3, 1), (3, 3), (4, 0), (4, 3)))]
assert placement2.point_multiplex_obs(gc1, 0) == [GriddedChord(Chord((0, 1, 1, 2, 0, 2)), ((0, 0), (1, 1), (3, 1), (3, 3), (4, 0), (4, 3)))]
'''
assert placement1.added_reqs(gc1, 5, 3, False) == [[GriddedChord(Chord((0,)), ((3,2),))],
                                                   [GriddedChord(Chord((0,)), ((1,2),))], 
                                                   [GriddedChord(Chord((0,1,1,2,0,2)), ((0,0),(0,0),(1,0),(1,2),(2,0),(3,2)))]]
assert placement1.added_reqs(gc1, 3, 5, True) == [[GriddedChord(Chord((0,)), ((2,2),))],
                                                  [GriddedChord(Chord((0,)), ((4,2),))], 
                                                  [GriddedChord(Chord((0,1,1,2,0,2)), ((0,0),(0,0),(1,0),(2,2),(4,0),(4,2)))]]
assert placement1.added_reqs(gc1, 0, 4, True) == [[GriddedChord(Chord((0,)), ((1,1),))],
                                                  [GriddedChord(Chord((0,)), ((4,1),))], 
                                                  [GriddedChord(Chord((0,1,1,2,0,2)), ((1,1),(2,2),(3,2),(3,3),(4,1),(4,3)))]]

assert placement1.point_dir_obs((2, 1), DIR_NORTH) == [GriddedChord(Chord((0,)), ((0, 3),)),
                                                       GriddedChord(Chord((0,)), ((1, 3),)),
                                                       GriddedChord(Chord((0,)), ((2, 3),)),
                                                       GriddedChord(Chord((0,)), ((3, 3),)),
                                                       GriddedChord(Chord((0,)), ((4, 3),))]
assert set(placement1ownrow.point_dir_obs((2, 0), DIR_SOUTH)) == set([GriddedChord(Chord((0,)), ((0, 0),)),
                                                       GriddedChord(Chord((0,)), ((1, 0),)),
                                                       GriddedChord(Chord((0,)), ((2, 0),))])
assert placement1.point_dir_obs((2, 0), DIR_EAST) == [GriddedChord(Chord((0,)), ((4, 0),)),
                                                       GriddedChord(Chord((0,)), ((4, 1),)),
                                                       GriddedChord(Chord((0,)), ((4, 2),)),
                                                       GriddedChord(Chord((0,)), ((4, 3),)),]
assert placement1owncol.point_dir_obs((0, 0), DIR_WEST) == [GriddedChord(Chord((0,)), ((0, 0),)),
                                                       GriddedChord(Chord((0,)), ((0, 1),)),]

print(placement_nc.place_point(GriddedChord(Chord((0,0)), ((0,0), (0,0))), 3, True))
print()
place = placement_nc.place_chord(GriddedChord(Chord((0,0)), ((0,0), (0,0))), 3)
place._remove_empty_rows_and_cols()
print(place)

