import _direct_run_bootstrap

from src.algorithms.factor import Factor
from src.common.chords import Chord, GriddedChord
from src.common.tiling import Tiling
from src.strategies.requirement_insertion import RequirementInsertionStrategy

gc_single_00_00 = GriddedChord(Chord((0, 0)), ((0, 0), (0, 0)))
gc_single_00_10 = GriddedChord(Chord((0, 0)), ((0, 0), (1, 0)))
gc_sc_disjoint = GriddedChord.single_cell(Chord((0, 0, 1, 1)), (0, 0))
gc_crossed = GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 1), (1, 0), (1, 1)))


def test_linkage_only_cells_are_active():
    tiling = Tiling((), (), (((0, 0), (1, 0)),), derive_empty=True)
    assert (0, 0) in tiling.active_cells
    assert (1, 0) in tiling.active_cells


def test_linkage_connectivity_in_contains():
    t_lk_single = Tiling((), (), (((0, 0),),), derive_empty=False, simplify=False)
    assert t_lk_single.contains(gc_single_00_00)
    assert not t_lk_single.contains(gc_single_00_10)
    assert not t_lk_single.contains(gc_sc_disjoint)


def test_hash_includes_linkages():
    t1 = Tiling((), (), (((0, 0),),), derive_empty=False)
    t2 = Tiling((), (), (((0, 0), (1, 0)),), derive_empty=False)
    assert hash(t1) != hash(t2)


def test_add_linkage():
    base = Tiling((), (), (), derive_empty=False)
    linked = base.add_linkage(((0, 0), (1, 0)))
    assert linked.linkages == (((0, 0), (1, 0)),)


def test_requirement_insertion_preserves_linkages():
    linkages = (((0, 0), (1, 0)),)
    parent = Tiling((), (), linkages, derive_empty=False)
    req = (GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),)
    strat = RequirementInsertionStrategy(req)
    avoid_child, contain_child = strat.decomposition_function(parent)
    assert avoid_child.linkages == linkages
    assert contain_child.linkages == linkages


def test_factor_unites_linkage_cells():
    tiling = Tiling(
        (),
        (),
        (((0, 0), (1, 0)),),
        derive_empty=False,
    )
    factor = Factor(tiling)
    assert not factor.factorable()


def test_factor_keeps_linkage_in_child():
    tiling = Tiling(
        (),
        (),
        (((0, 0), (1, 0)),),
        derive_empty=False,
        simplify=False,
    )
    factors = Factor(tiling).factors()
    assert len(factors) == 1
    assert factors[0].linkages == (((0, 0), (1, 0)),)


def test_sub_tiling_filters_partial_linkages():
    tiling = Tiling(
        (),
        (),
        (((0, 0), (1, 0)),),
        derive_empty=False,
    )
    child = tiling.sub_tiling(((0, 0),))
    assert child.linkages == ()


def test_json_roundtrip_with_linkages():
    tiling = Tiling((), ((gc_single_00_00,),), (((0, 0), (0, 1)),), derive_empty=False)
    assert Tiling.from_dict(tiling.to_jsonable()) == tiling


def test_simplify_deletes_singleton_linkage():
    tiling = Tiling((), (), (((0, 0),),), derive_empty=False, simplify=True)
    assert tiling.linkages == ()


def test_simplify_deletes_subset_linkage():
    tiling = Tiling(
        (),
        (),
        (((0, 0),), ((0, 0), (1, 0))),
        derive_empty=False,
        simplify=True,
    )
    assert tiling.linkages == (((0, 0), (1, 0)),)


def test_simplify_deletes_union_linkage():
    tiling = Tiling(
        (),
        (),
        (
            ((0, 0), (1, 0)),
            ((1, 0), (2, 0)),
            ((0, 0), (1, 0), (2, 0)),
        ),
        derive_empty=False,
        simplify=True,
        expand=False,
    )
    assert tiling.linkages == (((0, 0), (1, 0)), ((1, 0), (2, 0)))
