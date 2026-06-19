import _direct_run_bootstrap 

from src.algorithms.simplify import SimplifyObstructionsAndRequirements

from src.common.chords import GriddedChord, Chord
from src.common.tiling import Tiling

ob_containing_ob = SimplifyObstructionsAndRequirements(
    (GriddedChord(Chord((0, 1, 0, 2, 1, 2)), ((0, 0),) * 6), GriddedChord(Chord((0, 1, 0, 1)), ((0, 0),) * 4)),
    (),
    (1, 1),
    ((0, 0),),
    (),
)
ob_containing_ob.remove_redundant_obstructions()
assert ob_containing_ob.obstructions == (GriddedChord(Chord((0, 1, 0, 1)), ((0,0),)*4),)

req_containing_ob = SimplifyObstructionsAndRequirements(
    (GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),),
    ((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),),),
    (1, 1),
    ((0, 0),),
    (),
)
req_containing_ob.simplify()
assert req_containing_ob.obstructions == (GriddedChord(Chord((0, 0)), ((0,0), (0,0))),)
assert req_containing_ob.requirements == ((),)

all_from_21 = SimplifyObstructionsAndRequirements(
    (GriddedChord(Chord((0,)), ((1, 0),)),
     GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
     GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))),
     GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (2, 0))),
     GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
     GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (2, 0), (2, 0), (2, 0))),
     GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (2, 0))),
     GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
     GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (2, 0), (2, 0), (2, 0))),
     GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0))),
     GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (2, 0), (2, 0), (2, 0))),
     GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))),
     GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))),
     GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))),
     GriddedChord(Chord((1, 0)), ((2, 0), (2, 0)))),
    ((GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),),),
    (3, 1),
    ((0, 0), (1, 0), (2, 0)),
    ((1, 0),),
)

all_from_21.simplify()

req_containing_req = SimplifyObstructionsAndRequirements(
    (),
    ((GriddedChord(Chord((0, 1, 0, 2, 1, 2)), ((0, 0),) * 6),
      GriddedChord(Chord((0, 1, 0, 1)), ((0, 0),) * 4),),),
    (1, 1),
    ((0, 0),),
    (),
)
req_containing_req.remove_redundant_requirements()
assert req_containing_req.requirements == ((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (0, 0), (0, 0))),),)

reqlist_containing_reqlist = SimplifyObstructionsAndRequirements(
    (),
    ((GriddedChord(Chord((0, 1, 0, 2, 1, 2)), ((0, 0),) * 6),),
     (GriddedChord(Chord((0, 1, 0, 1)), ((0, 0),) * 4),),),
    (1, 1),
    ((0, 0),),
    (),
)
reqlist_containing_reqlist.remove_redundant_lists_requirements()
assert reqlist_containing_reqlist.requirements == (((GriddedChord(Chord((0, 1, 0, 2, 1, 2)), ((0, 0), )*6),),))

reqlist_containing_ob = SimplifyObstructionsAndRequirements(
    (GriddedChord(Chord((0, 1, 0, 1)), ((0, 0),) * 4),),
    ((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0),) * 4),
      GriddedChord(Chord((0, 1, 1, 0)), ((0, 0),) * 4),),),
    (1, 1),
    ((0, 0),),
    (),
)
reqlist_containing_ob.remove_redundant_requirements()
assert reqlist_containing_ob.requirements == ((GriddedChord(Chord((0,1,1,0)), ((0,0),)*4),),)
assert reqlist_containing_ob.obstructions == (GriddedChord(Chord((0, 1, 0, 1)), ((0,0),)*4),)

# Linkage deletion via simplify (Nabergall §2.2.1)
_row_3 = ((0, 0), (1, 0), (2, 0))
_dims_3x1 = (3, 1)
_empty = ()


def test_simplify_removes_singleton_linkage():
    algo = SimplifyObstructionsAndRequirements(
        (),
        (),
        (1, 1),
        ((0, 0),),
        _empty,
        linkages=(((0, 0),),),
    )
    algo.simplify()
    assert algo.linkages == ()


def test_simplify_removes_subset_linkage():
    algo = SimplifyObstructionsAndRequirements(
        (),
        (),
        (2, 1),
        ((0, 0), (1, 0)),
        _empty,
        linkages=(((0, 0),), ((0, 0), (1, 0))),
    )
    algo.simplify()
    assert algo.linkages == (((0, 0), (1, 0)),)


def test_simplify_removes_union_linkage():
    algo = SimplifyObstructionsAndRequirements(
        (),
        (),
        _dims_3x1,
        _row_3,
        _empty,
        linkages=(
            ((0, 0), (1, 0)),
            ((1, 0), (2, 0)),
            ((0, 0), (1, 0), (2, 0)),
        ),
    )
    algo.simplify()
    assert algo.linkages == (((0, 0), (1, 0)), ((1, 0), (2, 0)))

