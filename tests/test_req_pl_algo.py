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
            GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),
            GriddedChord(Chord((0, 0)), ((1, 0), (2, 0))),
            GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))),
            GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))),
            GriddedChord(Chord((1, 0)), ((2, 0), (2, 0))),
            GriddedChord(Chord((1, 0)), ((0, 1), (1, 1))),
        ],
        requirements=[
            [GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))],
            [GriddedChord(Chord((0, 0)), ((1, 1), (2, 1))),
             GriddedChord(Chord((0, 0)), ((2, 1), (2, 1))),
            ],
        ],
    )

tiling2 = Tiling(
        obstructions=[
            GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), ) * 4),
        ],
        requirements=[
            [GriddedChord(Chord((0, 0)), ((0, 0), (2, 0)))],
            [GriddedChord(Chord((0, 1)), ((0, 0), (1, 1)))],
        ],
    )

placement1 = RequirementPlacement(tiling1)
placement1ownrow = RequirementPlacement(tiling1, own_row=True, own_col=False)
placement1owncol = RequirementPlacement(tiling1, own_row=False, own_col=True)

placement2 = RequirementPlacement(tiling2)
placement2ownrow = RequirementPlacement(tiling2, own_row=True, own_col=False)
placement2owncol = RequirementPlacement(tiling2, own_row=False, own_col=True)

placement_nc = RequirementPlacement(non_crossing)
placement_nc_ownroq = RequirementPlacement(tiling2, own_row=True, own_col=False)
placement_nc_owncol = RequirementPlacement(tiling2, own_row=False, own_col=True)

t = Tiling(
        obstructions=[GriddedChord(Chord((0, 1, 1, 0)), ((1, 1),) * 4)],
        requirements=[[GriddedChord(Chord((0,)), ((0, 0),))]],
    )
placement_only_west = RequirementPlacement(t, dirs=[DIR_WEST])

gc1 = GriddedChord(Chord((0,1,1,2,0,2)), ((0, 0), (0, 0), (1, 0), (1, 1), (1, 0), (1,1)))

print(tiling1.point_cells)

#print(placement_nc._tiling)
#tilings = placement_nc.col_placement(1, DIR_WEST)
#print(len(tilings))
#for t in tilings:
#    print(t)

print(placement_only_west._tiling_point_col_cells())
assert placement2._tiling_point_col_cells() ==  frozenset([])
assert placement1._tiling_point_col_cells() == frozenset([])

t_0 = Tiling((GriddedChord(Chord((0,1)), ((0,1), (0,1))), GriddedChord(Chord((0,0)), ((0,1), (0,1))), GriddedChord(Chord((1,0)), ((0,1), (0,1))), GriddedChord(Chord((0,1)), ((0,0), (0,0))), GriddedChord(Chord((0,0)), ((0,0), (0,0))), GriddedChord(Chord((1,0)), ((0,0), (0,0)))),
           ((GriddedChord(Chord((0,)), ((0,0),)),),(GriddedChord(Chord((0,)), ((0,1),)),)))

req_pl = RequirementPlacement(t_0)

print(req_pl._tiling_point_col_cells())
print(req_pl._tiling_point_row_cells())

'''
# ------------------------------------------------------------
#       Tests for RequirementPlacement Class
# ------------------------------------------------------------

def test_col_placement(placement1: "RequirementPlacement"):
    print(placement1._tiling)
    tilings = placement1.col_placement(1, DIR_WEST)
    assert len(tilings) == 2
    assert all(isinstance(t, Tiling) for t in tilings)


def test_row_placement(placement1: "RequirementPlacement"):
    print(placement1._tiling)
    tilings = placement1.row_placement(1, DIR_NORTH)
    assert len(tilings) == 1
    assert all(isinstance(t, Tiling) for t in tilings)


def test_empty_row(placement1):
    t = Tiling(
        obstructions=(
            GriddedChord((0,), ((0, 1),)),
            GriddedChord((0,), ((1, 0),)),
            GriddedChord((0,), ((2, 0),)),
            GriddedChord((0,), ((3, 1),)),
            GriddedChord((1, 0), ((2, 1), (2, 1))),
            GriddedChord((0, 1, 2), ((1, 1), (1, 1), (1, 1))),
            GriddedChord((2, 0, 1), ((3, 0), (3, 0), (3, 0))),
            GriddedChord((2, 1, 0), ((0, 0), (0, 0), (0, 0))),
        ),
    )
    assert placement1.empty_row(1) == t


def test_empty_col(placement1):
    t = Tiling(
        obstructions=(
            GriddedChord((0,), ((0, 1),)),
            GriddedChord((0,), ((0, 2),)),
            GriddedChord((0,), ((1, 0),)),
            GriddedChord((0,), ((1, 1),)),
            GriddedChord((0,), ((2, 1),)),
            GriddedChord((0,), ((2, 2),)),
            GriddedChord((1, 0), ((1, 2), (1, 2))),
            GriddedChord((2, 0, 1), ((2, 0), (2, 0), (2, 0))),
            GriddedChord((2, 1, 0), ((0, 0), (0, 0), (0, 0))),
        ),
        requirements=(),
    )
    print(t)
    print(placement1._tiling)
    assert placement1.empty_col(1) == t


def test_point_translation(gp1, placement1: RequirementPlacement, placement1owncol, placement1ownrow):
    assert placement1._point_translation(gp1, 2, (0, 3)) == (3, 1)
    assert placement1._point_translation(gp1, 2, (1, 2)) == (3, 3)
    assert placement1._point_translation(gp1, 2, (2, 2)) == (3, 3)
    assert placement1._point_translation(gp1, 2, (3, 0)) == (1, 3)
    assert placement1._point_translation(gp1, 2, (4, 4)) == (1, 1)

    assert placement1owncol._point_translation(gp1, 2, (0, 3)) == (3, 1)
    assert placement1owncol._point_translation(gp1, 2, (1, 2)) == (3, 1)
    assert placement1owncol._point_translation(gp1, 2, (2, 2)) == (3, 1)
    assert placement1owncol._point_translation(gp1, 2, (3, 0)) == (1, 1)
    assert placement1owncol._point_translation(gp1, 2, (4, 4)) == (1, 1)

    assert placement1ownrow._point_translation(gp1, 2, (0, 3)) == (1, 1)
    assert placement1ownrow._point_translation(gp1, 2, (1, 2)) == (1, 3)
    assert placement1ownrow._point_translation(gp1, 2, (2, 2)) == (1, 3)
    assert placement1ownrow._point_translation(gp1, 2, (3, 0)) == (1, 3)
    assert placement1ownrow._point_translation(gp1, 2, (4, 4)) == (1, 1)


def test_gridded_perm_translation(gp1, placement1, placement1owncol, placement1ownrow):
    assert placement1._gridded_perm_translation(gp1, (0, 3)) == GriddedChord(
        (3, 1, 2, 0, 4), ((2, 3), (2, 0), (3, 1), (3, 0), (3, 3))
    )
    assert placement1._gridded_perm_translation(gp1, (1, 1)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (2, 2), (3, 3), (3, 0), (3, 3))
    )
    assert placement1._gridded_perm_translation(gp1, (2, 2)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 0), (3, 3), (3, 0), (3, 3))
    )
    assert placement1._gridded_perm_translation(gp1, (3, 0)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 2), (1, 3), (3, 2), (3, 3))
    )
    assert placement1._gridded_perm_translation(gp1, (4, 4)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (3, 3))
    )
    assert placement1owncol._gridded_perm_translation(gp1, (0, 3)) == GriddedChord(
        (3, 1, 2, 0, 4), ((2, 1), (2, 0), (3, 1), (3, 0), (3, 1))
    )
    assert placement1owncol._gridded_perm_translation(gp1, (1, 1)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (2, 0), (3, 1), (3, 0), (3, 1))
    )
    assert placement1owncol._gridded_perm_translation(gp1, (2, 2)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (3, 1), (3, 0), (3, 1))
    )
    assert placement1owncol._gridded_perm_translation(gp1, (3, 0)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (3, 0), (3, 1))
    )
    assert placement1owncol._gridded_perm_translation(gp1, (4, 4)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (3, 1))
    )
    assert placement1ownrow._gridded_perm_translation(gp1, (0, 3)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 1), (1, 0), (1, 3))
    )
    assert placement1ownrow._gridded_perm_translation(gp1, (1, 1)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 2), (1, 3), (1, 0), (1, 3))
    )
    assert placement1ownrow._gridded_perm_translation(gp1, (2, 2)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 3), (1, 0), (1, 3))
    )
    assert placement1ownrow._gridded_perm_translation(gp1, (3, 0)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 2), (1, 3), (1, 2), (1, 3))
    )
    assert placement1ownrow._gridded_perm_translation(gp1, (4, 4)) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (1, 3))
    )


def test_gridded_perm_translation_with_point(
    gp1, placement1, placement1owncol, placement1ownrow
):
    assert placement1._gridded_perm_translation_with_point(gp1, 0) == GriddedChord(
        (3, 1, 2, 0, 4), ((1, 2), (2, 0), (3, 1), (3, 0), (3, 3))
    )
    assert placement1._gridded_perm_translation_with_point(gp1, 1) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (1, 1), (3, 3), (3, 0), (3, 3))
    )
    assert placement1._gridded_perm_translation_with_point(gp1, 2) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 0), (2, 2), (3, 0), (3, 3))
    )
    assert placement1._gridded_perm_translation_with_point(gp1, 3) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 2), (1, 3), (2, 1), (3, 3))
    )
    assert placement1._gridded_perm_translation_with_point(gp1, 4) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (2, 2))
    )
    assert placement1ownrow._gridded_perm_translation_with_point(gp1, 0) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 2), (0, 0), (1, 1), (1, 0), (1, 3))
    )
    assert placement1ownrow._gridded_perm_translation_with_point(gp1, 1) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 1), (1, 3), (1, 0), (1, 3))
    )
    assert placement1ownrow._gridded_perm_translation_with_point(gp1, 2) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 2), (1, 0), (1, 3))
    )
    assert placement1ownrow._gridded_perm_translation_with_point(gp1, 3) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 3), (0, 2), (1, 3), (1, 1), (1, 3))
    )
    assert placement1ownrow._gridded_perm_translation_with_point(gp1, 4) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (1, 2))
    )
    assert placement1owncol._gridded_perm_translation_with_point(gp1, 0) == GriddedChord(
        (3, 1, 2, 0, 4), ((1, 1), (2, 0), (3, 1), (3, 0), (3, 1))
    )
    assert placement1owncol._gridded_perm_translation_with_point(gp1, 1) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (1, 0), (3, 1), (3, 0), (3, 1))
    )
    assert placement1owncol._gridded_perm_translation_with_point(gp1, 2) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (2, 1), (3, 0), (3, 1))
    )
    assert placement1owncol._gridded_perm_translation_with_point(gp1, 3) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (2, 0), (3, 1))
    )
    assert placement1owncol._gridded_perm_translation_with_point(gp1, 4) == GriddedChord(
        (3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (2, 1))
    )


def test_placed_cell(placement1, placement1owncol, placement1ownrow):
    assert placement1._placed_cell((0, 0)) == (1, 1)
    assert placement1._placed_cell((3, 2)) == (4, 3)
    assert placement1owncol._placed_cell((9, 7)) == (10, 7)
    assert placement1owncol._placed_cell((2, 1)) == (3, 1)
    assert placement1ownrow._placed_cell((0, 4)) == (0, 5)
    assert placement1ownrow._placed_cell((4, 2)) == (4, 3)


def test_point_obstructions(placement1, placement1owncol, placement1ownrow):
    assert placement1._point_obstructions((0, 0)) == [
        GriddedChord((0, 1), ((1, 1), (1, 1))),
        GriddedChord((1, 0), ((1, 1), (1, 1))),
    ]
    assert placement1owncol._point_obstructions((0, 0)) == [
        GriddedChord((0, 1), ((1, 0), (1, 0))),
        GriddedChord((1, 0), ((1, 0), (1, 0))),
    ]
    assert placement1ownrow._point_obstructions((0, 0)) == [
        GriddedChord((0, 1), ((0, 1), (0, 1))),
        GriddedChord((1, 0), ((0, 1), (0, 1))),
    ]


def test_point_requirements(placement1, placement1owncol, placement1ownrow):
    assert placement1._point_requirements((2, 3)) == [[GriddedChord((0,), ((3, 4),))]]
    assert placement1ownrow._point_requirements((2, 3)) == [
        [GriddedChord((0,), ((2, 4),))]
    ]
    assert placement1owncol._point_requirements((2, 3)) == [
        [GriddedChord((0,), ((3, 3),))]
    ]


def test_stretch_gridded_perm(gp1, placement1, placement1owncol, placement1ownrow):
    assert set(placement1._stretch_gridded_perm(gp1, (0, 0))) == set(
        [
            GriddedChord((3, 1, 2, 0, 4), ((2, 3), (2, 2), (3, 3), (3, 2), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((2, 3), (2, 2), (3, 3), (3, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((2, 3), (2, 0), (3, 3), (3, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (2, 2), (3, 3), (3, 2), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (2, 2), (3, 3), (3, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (2, 0), (3, 3), (3, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 2), (3, 3), (3, 2), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 2), (3, 3), (3, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (3, 3), (3, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (1, 1), (3, 3), (3, 0), (3, 3))),
        ]
    )
    assert set(placement1owncol._stretch_gridded_perm(gp1, (1, 0))) == set(
        [
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (0, 0), (3, 1), (3, 0), (3, 1))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (3, 0), (3, 1))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (3, 1))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (1, 1))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (2, 0), (3, 1))),
        ]
    )
    assert set(placement1ownrow._stretch_gridded_perm(gp1, (1, 1))) == set(
        [
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 3), (1, 0), (1, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 1), (1, 0), (1, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (1, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (1, 1))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 2), (1, 0), (1, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (0, 0), (1, 1), (1, 0), (1, 2))),
        ]
    )


def test_stretch_gridded_perms(placement1, placement1owncol, placement1ownrow):
    gps = [
        GriddedChord((0, 1), [(0, 0), (1, 1)]),
        GriddedChord((0, 1), [(1, 1), (2, 2)]),
    ]
    for p in (placement1, placement1ownrow, placement1owncol):
        assert set(p._stretch_gridded_perms(gps, (1, 1))) == set(
            chain.from_iterable(p._stretch_gridded_perm(gp, (1, 1)) for gp in gps)
        )


def test_stretched_obstructions(placement1, placement1owncol, placement1ownrow):
    orig_obs = placement1._tiling.obstructions
    assert sorted(placement1.stretched_obstructions((1, 1))) == sorted(
        placement1._stretch_gridded_perms(orig_obs, (1, 1))
    )
    assert sorted(placement1owncol.stretched_obstructions((1, 1))) == sorted(
        placement1owncol._stretch_gridded_perms(orig_obs, (1, 1))
    )
    assert sorted(placement1ownrow.stretched_obstructions((1, 1))) == sorted(
        placement1ownrow._stretch_gridded_perms(orig_obs, (1, 1))
    )


def test_stretched_requirements(placement1, placement1owncol, placement1ownrow):
    orig_reqs = placement1._tiling.requirements
    assert sorted(placement1.stretched_requirements((1, 1))) == sorted(
        placement1._stretch_gridded_perms(orig_reqs, (1, 1))
    )
    orig_reqs = placement1owncol._tiling.requirements
    assert sorted(placement1owncol.stretched_requirements((1, 1))) == sorted(
        placement1owncol._stretch_gridded_perms(orig_reqs, (1, 1))
    )
    orig_reqs = placement1ownrow._tiling.requirements
    assert sorted(placement1ownrow.stretched_requirements((1, 1))) == sorted(
        placement1ownrow._stretch_gridded_perms(orig_reqs, (1, 1))
    )


def test_stretched_obstructions_and_assumptions(
    placement1, placement1owncol, placement1ownrow
):
    obs, reqs, _ = placement1._stretched_obstructions_requirements_and_assumptions(
        (1, 1)
    )
    assert set(obs) == set(
        placement1.stretched_obstructions((1, 1))
        + [
            GriddedChord.single_cell((0, 1), (2, 2)),
            GriddedChord.single_cell((1, 0), (2, 2)),
        ]
    )
    assert sorted(reqs) == sorted(
        placement1.stretched_requirements((1, 1)) + [[GriddedChord((0,), ((2, 2),))]]
    )
    (
        obs,
        reqs,
        _,
    ) = placement1ownrow._stretched_obstructions_requirements_and_assumptions((1, 1))
    assert set(obs) == set(
        placement1ownrow.stretched_obstructions((1, 1))
        + [
            GriddedChord.single_cell((0, 1), (1, 2)),
            GriddedChord.single_cell((1, 0), (1, 2)),
        ]
    )
    assert sorted(reqs) == sorted(
        placement1ownrow.stretched_requirements((1, 1))
        + [[GriddedChord((0,), ((1, 2),))]]
    )
    (
        obs,
        reqs,
        _,
    ) = placement1owncol._stretched_obstructions_requirements_and_assumptions((1, 1))
    assert set(obs) == set(
        placement1owncol.stretched_obstructions((1, 1))
        + [
            GriddedChord.single_cell((0, 1), (2, 1)),
            GriddedChord.single_cell((1, 0), (2, 1)),
        ]
    )
    assert sorted(reqs) == sorted(
        placement1owncol.stretched_requirements((1, 1))
        + [[GriddedChord((0,), ((2, 1),))]]
    )


def farther(placement1):
    assert placement1._farther((0, 0), (2, 0), DIR_EAST) is False
    assert placement1._farther((0, 0), (2, 0), DIR_NORTH) is False
    assert placement1._farther((0, 0), (2, 0), DIR_WEST) is True
    assert placement1._farther((0, 0), (2, 0), DIR_SOUTH) is False

    assert placement1._farther((2, 3), (2, 0), DIR_EAST) is False
    assert placement1._farther((2, 3), (2, 0), DIR_NORTH) is True
    assert placement1._farther((2, 3), (2, 0), DIR_WEST) is False
    assert placement1._farther((2, 3), (2, 0), DIR_SOUTH) is False

    assert placement1._farther((1, 1), (3, 4), DIR_EAST) is False
    assert placement1._farther((1, 1), (3, 4), DIR_NORTH) is False
    assert placement1._farther((1, 1), (3, 4), DIR_WEST) is True
    assert placement1._farther((1, 1), (3, 4), DIR_SOUTH) is True

    assert placement1._farther((1, 5), (3, 4), DIR_EAST) is False
    assert placement1._farther((1, 5), (3, 4), DIR_NORTH) is True
    assert placement1._farther((1, 5), (3, 4), DIR_WEST) is True
    assert placement1._farther((1, 5), (3, 4), DIR_SOUTH) is False

    assert placement1._farther((2, 2), (1, 1), DIR_EAST) is True
    assert placement1._farther((2, 2), (1, 1), DIR_NORTH) is True
    assert placement1._farther((2, 2), (1, 1), DIR_WEST) is False
    assert placement1._farther((2, 2), (1, 1), DIR_SOUTH) is False


def test_forced_obstructions_from_patt(
    gp1, placement1, placement1owncol, placement1ownrow
):
    assert set(
        placement1.forced_obstructions_from_requirement(
            (gp1,), (2,), gp1.pos[2], DIR_NORTH
        )
    ) == set(
        [
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (3, 3), (3, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 3), (3, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 3), (1, 0), (3, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 3), (1, 0), (1, 3))),
        ]
    )

    assert set(
        placement1owncol.forced_obstructions_from_requirement(
            (gp1,), (1,), gp1.pos[1], DIR_EAST
        )
    ) == set(
        [
            GriddedChord((3, 1, 2, 0, 4), ((2, 1), (2, 0), (3, 1), (3, 0), (3, 1))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 1), (2, 0), (3, 1), (3, 0), (3, 1))),
        ]
    )

    assert set(
        placement1ownrow.forced_obstructions_from_requirement(
            (gp1,), (3,), gp1.pos[3], DIR_SOUTH
        )
    ) == set(
        [
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 2), (1, 3), (1, 0), (1, 3))),
            GriddedChord((3, 1, 2, 0, 4), ((0, 3), (0, 0), (1, 3), (1, 0), (1, 3))),
        ]
    )


def test_forced_obstructions_from_list(
    gp1, placement1, placement1owncol, placement1ownrow
):
    req_list_row = [
        GriddedChord((0,), ((0, 0),)),
        GriddedChord((0,), ((1, 0),)),
    ]
    assert set(
        placement1.forced_obstructions_from_requirement(
            req_list_row, (0, 0), (0, 0), DIR_NORTH
        )
    ) == set(
        [
            GriddedChord((0,), ((0, 2),)),
            GriddedChord((0,), ((2, 2),)),
            GriddedChord((0,), ((3, 2),)),
        ]
    )
    assert set(
        placement1.forced_obstructions_from_requirement(
            req_list_row, (0, 0), (0, 0), DIR_SOUTH
        )
    ) == set(
        [
            GriddedChord((0,), ((0, 0),)),
            GriddedChord((0,), ((2, 0),)),
            GriddedChord((0,), ((3, 0),)),
        ]
    )
    assert set(
        placement1.forced_obstructions_from_requirement(
            req_list_row, (0, 0), (1, 0), DIR_NORTH
        )
    ) == set(
        [
            GriddedChord((0,), ((0, 2),)),
            GriddedChord((0,), ((1, 2),)),
            GriddedChord((0,), ((3, 2),)),
        ]
    )
    assert set(
        placement1.forced_obstructions_from_requirement(
            req_list_row, (0, 0), (1, 0), DIR_SOUTH
        )
    ) == set(
        [
            GriddedChord((0,), ((0, 0),)),
            GriddedChord((0,), ((1, 0),)),
            GriddedChord((0,), ((3, 0),)),
        ]
    )
    assert set(
        placement1ownrow.forced_obstructions_from_requirement(
            req_list_row, (0, 0), (0, 0), DIR_NORTH
        )
    ) == set([GriddedChord((0,), ((0, 2),)), GriddedChord((0,), ((1, 2),))])
    assert set(
        placement1ownrow.forced_obstructions_from_requirement(
            req_list_row, (0, 0), (0, 0), DIR_SOUTH
        )
    ) == set([GriddedChord((0,), ((0, 0),)), GriddedChord((0,), ((1, 0),))])
    assert set(
        placement1ownrow.forced_obstructions_from_requirement(
            req_list_row, (0, 0), (1, 0), DIR_NORTH
        )
    ) == set([GriddedChord((0,), ((0, 2),)), GriddedChord((0,), ((1, 2),))])
    assert set(
        placement1ownrow.forced_obstructions_from_requirement(
            req_list_row, (0, 0), (1, 0), DIR_SOUTH
        )
    ) == set([GriddedChord((0,), ((0, 0),)), GriddedChord((0,), ((1, 0),))])

    req_list_col = [
        GriddedChord((0,), ((0, 0),)),
        GriddedChord((0,), ((0, 1),)),
    ]
    assert set(
        placement1.forced_obstructions_from_requirement(
            req_list_col, (0, 0), (0, 0), DIR_EAST
        )
    ) == set(
        [
            GriddedChord((0,), ((2, 0),)),
            GriddedChord((0,), ((2, 2),)),
            GriddedChord((0,), ((2, 3),)),
        ]
    )
    assert set(
        placement1.forced_obstructions_from_requirement(
            req_list_col, (0, 0), (0, 0), DIR_WEST
        )
    ) == set(
        [
            GriddedChord((0,), ((0, 0),)),
            GriddedChord((0,), ((0, 2),)),
            GriddedChord((0,), ((0, 3),)),
        ]
    )
    assert set(
        placement1.forced_obstructions_from_requirement(
            req_list_col, (0, 0), (0, 1), DIR_EAST
        )
    ) == set(
        [
            GriddedChord((0,), ((2, 0),)),
            GriddedChord((0,), ((2, 1),)),
            GriddedChord((0,), ((2, 3),)),
        ]
    )
    assert set(
        placement1.forced_obstructions_from_requirement(
            req_list_col, (0, 0), (0, 1), DIR_WEST
        )
    ) == set(
        [
            GriddedChord((0,), ((0, 0),)),
            GriddedChord((0,), ((0, 1),)),
            GriddedChord((0,), ((0, 3),)),
        ]
    )
    assert set(
        placement1owncol.forced_obstructions_from_requirement(
            req_list_col, (0, 0), (0, 0), DIR_EAST
        )
    ) == set([GriddedChord((0,), ((2, 0),)), GriddedChord((0,), ((2, 1),))])
    assert set(
        placement1owncol.forced_obstructions_from_requirement(
            req_list_col, (0, 0), (0, 0), DIR_WEST
        )
    ) == set([GriddedChord((0,), ((0, 0),)), GriddedChord((0,), ((0, 1),))])
    assert set(
        placement1owncol.forced_obstructions_from_requirement(
            req_list_col, (0, 0), (0, 1), DIR_EAST
        )
    ) == set([GriddedChord((0,), ((2, 0),)), GriddedChord((0,), ((2, 1),))])
    assert set(
        placement1owncol.forced_obstructions_from_requirement(
            req_list_col, (0, 0), (0, 1), DIR_WEST
        )
    ) == set([GriddedChord((0,), ((0, 0),)), GriddedChord((0,), ((0, 1),))])
'''