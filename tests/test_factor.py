import _direct_run_bootstrap

from comb_spec_searcher.strategies.strategy import strategy_from_dict

from src.algorithms.factor import Factor
from src.common.assumptions import TrackingAssumption
from src.common.chords import Chord, GriddedChord
from src.common.tiling import Tiling
from src.strategies.factor import FactorFactory, FactorStrategy


def crossed_chord(pos):
    return GriddedChord(Chord((0, 1, 0, 1)), (pos,) * 4)


def single_chord(pos_left, pos_right):
    return GriddedChord(Chord((0, 0)), (pos_left, pos_right))


def _make_non_crossing() -> Tiling:
    return Tiling(
        obstructions=sorted(
            (
                GriddedChord(Chord((0,)), ((0, 1),)),
                GriddedChord(Chord((0,)), ((0, 2),)),
                GriddedChord(Chord((0,)), ((1, 0),)),
                GriddedChord(Chord((0,)), ((1, 1),)),
                GriddedChord(Chord((0,)), ((2, 1),)),
                GriddedChord(Chord((0,)), ((2, 2),)),
                GriddedChord(Chord((0,)), ((3, 0),)),
                GriddedChord(Chord((0,)), ((3, 2),)),
                GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))),
                GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))),
                GriddedChord(Chord((0, 1)), ((0, 0), (2, 0))),
                GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))),
                GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))),
                GriddedChord(Chord((1, 0)), ((0, 0), (2, 0))),
                GriddedChord(Chord((1, 0)), ((2, 0), (2, 0))),
                GriddedChord(Chord((0, 1, 0, 1)), ((1, 2),) * 4),
                GriddedChord(Chord((0, 1, 0, 1)), ((3, 1),) * 4),
            )
        ),
        requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),),),
    )


def _disjoint_two_cell_tiling() -> Tiling:
    """Two active cells on distinct rows and columns — should factor."""
    return Tiling(
        obstructions=(
            crossed_chord((0, 0)),
            crossed_chord((1, 1)),
        ),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )


def _disjoint_three_cell_tiling() -> Tiling:
    return Tiling(
        obstructions=(
            crossed_chord((0, 0)),
            crossed_chord((1, 1)),
            crossed_chord((2, 2)),
        ),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )


# ---------------------------------------------------------------------------
# algorithms/factor.py
# ---------------------------------------------------------------------------


def test_factor_decomposes_into_multiple_factors():
    tiling = _make_non_crossing()
    algo = Factor(tiling)
    assert algo.factorable()

    factors = algo.factors()
    assert isinstance(factors, tuple)
    assert len(factors) >= 2
    assert all(isinstance(t, Tiling) for t in factors)


def test_factor_not_factorable_when_single_component():
    tiling = Tiling(
        obstructions=(crossed_chord((0, 0)),),
        simplify=False,
        expand=False,
    )
    algo = Factor(tiling)
    assert not algo.factorable()
    assert len(algo.get_components()) == 1
    assert len(algo.factors()) == 1


def test_factor_components_are_cached():
    tiling = _disjoint_two_cell_tiling()
    algo = Factor(tiling)
    first = algo.get_components()
    second = algo.get_components()
    assert first is second
    assert len(first) == 2
    assert {frozenset(c) for c in first} == {frozenset([(0, 0)]), frozenset([(1, 1)])}


def test_factor_unites_same_row_or_column():
    same_row = Tiling(
        obstructions=(crossed_chord((0, 0)), crossed_chord((1, 0))),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    same_col = Tiling(
        obstructions=(crossed_chord((0, 0)), crossed_chord((0, 1))),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    assert not Factor(same_row).factorable()
    assert not Factor(same_col).factorable()


def test_factor_unites_cells_sharing_an_obstruction():
    tiling = Tiling(
        obstructions=(
            crossed_chord((0, 0)),
            crossed_chord((1, 1)),
            single_chord((0, 0), (1, 1)),
        ),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    assert not Factor(tiling).factorable()
    assert len(Factor(tiling).get_components()) == 1


def test_factor_unites_cells_sharing_a_requirement():
    tiling = Tiling(
        obstructions=(
            crossed_chord((0, 0)),
            crossed_chord((1, 1)),
        ),
        requirements=((single_chord((0, 0), (1, 1)),),),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    assert not Factor(tiling).factorable()


def test_factor_unites_linkage_cells():
    tiling = Tiling(
        (),
        (),
        (((0, 0), (1, 1)),),
        derive_empty=False,
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    # Linkage alone can activate cells; they must stay in one component.
    assert not Factor(tiling).factorable()
    factors = Factor(tiling).factors()
    assert len(factors) == 1
    assert factors[0].linkages == (((0, 0), (1, 1)),)


def test_factor_unites_cells_inside_tracked_assumption_chord():
    """Cells of one tracked multi-cell chord end up in the same factor."""
    ass = TrackingAssumption([single_chord((0, 0), (1, 1))])
    tiling = Tiling(
        obstructions=(
            crossed_chord((0, 0)),
            crossed_chord((1, 1)),
        ),
        assumptions=(ass,),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    assert not Factor(tiling).factorable()


def test_factor_factors_obs_and_reqs_split_correctly():
    tiling = Tiling(
        obstructions=(
            crossed_chord((0, 0)),
            crossed_chord((1, 1)),
        ),
        requirements=(
            (GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),),
            (GriddedChord(Chord((0, 0)), ((1, 1), (1, 1))),),
        ),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    algo = Factor(tiling)
    parts = algo._get_factors_obs_and_reqs()
    assert len(parts) == 2
    # Each factor keeps the obstruction/requirement supported on its cells.
    for obs, reqs, _assumps in parts:
        assert len(obs) == 1
        assert len(reqs) == 1
        cells = {c for ob in obs for c in ob.pos}
        assert all(req[0].pos[0] in cells for req in reqs)

    # Cached on second call.
    assert algo._get_factors_obs_and_reqs() is parts


def test_factor_empty_tiling_obs_and_reqs():
    tiling = Tiling(
        obstructions=(GriddedChord(Chord(()), ()),),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    assert tiling.is_empty()
    parts = Factor(tiling)._get_factors_obs_and_reqs()
    assert len(parts) == 1
    obs, reqs, assumps = parts[0]
    assert reqs == ()
    assert assumps == ()
    assert len(obs) == 1
    assert obs[0].patt == Chord(())


def test_factor_reducible_factorisations_groups_irreducible_factors():
    tiling = _disjoint_three_cell_tiling()
    algo = Factor(tiling)
    assert len(algo.get_components()) == 3

    reducible = list(algo.reducible_factorisations())
    # partitions_iterator yields non-trivial partitions of 3 parts: 3 of type 2+1
    assert len(reducible) == 3
    for factorization in reducible:
        assert len(factorization) == 2
        assert all(isinstance(t, Tiling) for t in factorization)

    # Two factors have no non-trivial reducible grouping.
    two = _disjoint_two_cell_tiling()
    assert list(Factor(two).reducible_factorisations()) == []


def test_factor_same_row_or_col_helper():
    assert Factor._same_row_or_col((0, 0), (1, 0))
    assert Factor._same_row_or_col((0, 0), (0, 1))
    assert not Factor._same_row_or_col((0, 0), (1, 1))


# ---------------------------------------------------------------------------
# strategies/factor.py
# ---------------------------------------------------------------------------


def test_factor_factory_yields_a_strategy():
    tiling = _make_non_crossing()
    factory = FactorFactory(True)
    strats = list(factory(tiling))
    assert strats
    assert isinstance(strats[0], FactorStrategy)


def test_factor_factory_smoke_on_small_tiling():
    tiling = Tiling(
        obstructions=(crossed_chord((0, 0)),),
        simplify=False,
        expand=False,
    )
    factory = FactorFactory(unions=False)
    assert list(factory(tiling)) == []


def test_factor_factory_yields_nothing_when_not_factorable():
    tiling = Tiling(
        obstructions=(crossed_chord((0, 0)), crossed_chord((1, 0))),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    assert list(FactorFactory()(tiling)) == []


def test_factor_factory_emits_irreducible_partition():
    tiling = _disjoint_two_cell_tiling()
    strats = list(FactorFactory(unions=False)(tiling))
    assert len(strats) == 1
    strat = strats[0]
    assert isinstance(strat, FactorStrategy)
    assert strat.workable is True
    assert strat.ignore_parent is True
    assert set(map(frozenset, strat.partition)) == {
        frozenset([(0, 0)]),
        frozenset([(1, 1)]),
    }

    children = strat.decomposition_function(tiling)
    assert len(children) == 2
    assert all(isinstance(c, Tiling) for c in children)
    assert all(len(c.active_cells) == 1 for c in children)


def test_factor_factory_unions_include_grouped_partitions():
    tiling = _disjoint_three_cell_tiling()
    without = list(FactorFactory(unions=False)(tiling))
    with_unions = list(FactorFactory(unions=True)(tiling))

    assert len(without) == 1
    # 3 union partitions (2+1) plus the fully split irreducible partition
    assert len(with_unions) == 4

    irreducible = without[0]
    assert irreducible in with_unions or irreducible.partition in {
        s.partition for s in with_unions
    }

    grouped = [s for s in with_unions if s.workable is False]
    assert len(grouped) == 3
    assert all(len(s.partition) == 2 for s in grouped)

    full = [s for s in with_unions if s.workable is True]
    assert len(full) == 1
    assert len(full[0].partition) == 3


def test_factor_strategy_formal_step_and_repr():
    tiling = _disjoint_two_cell_tiling()
    strat = next(FactorFactory()(tiling))
    assert "factor with partition" in strat.formal_step()
    assert "{(0, 0)}" in strat.formal_step()
    assert "{(1, 1)}" in strat.formal_step()
    assert str(strat) == "factor"
    assert "FactorStrategy" in repr(strat)
    assert "partition=" in repr(strat)


def test_factor_strategy_json_roundtrip():
    tiling = _disjoint_two_cell_tiling()
    strat = next(FactorFactory(ignore_parent=False, workable=False)(tiling))
    restored = strategy_from_dict(strat.to_jsonable())
    assert isinstance(restored, FactorStrategy)
    assert restored.partition == strat.partition
    assert restored.ignore_parent == strat.ignore_parent
    assert restored.workable == strat.workable


def test_factor_factory_json_roundtrip():
    factory = FactorFactory(unions=True, ignore_parent=False, workable=False)
    restored = strategy_from_dict(factory.to_jsonable())
    assert isinstance(restored, FactorFactory)
    assert restored.unions is True
    assert restored.ignore_parent is False
    assert restored.workable is False
    assert "unions of factor" in str(factory)
    assert "FactorFactory" in repr(factory)


def test_factor_strategy_forward_backward_map_roundtrip():
    tiling = _disjoint_two_cell_tiling()
    strat = next(FactorFactory()(tiling))
    children = strat.decomposition_function(tiling)

    parent_gc = GriddedChord(
        Chord((0, 0, 1, 1)),
        ((0, 0), (0, 0), (1, 1), (1, 1)),
    )
    assert tiling.contains(parent_gc)

    mapped = strat.forward_map(tiling, parent_gc, children)
    assert len(mapped) == 2
    assert all(isinstance(gc, GriddedChord) for gc in mapped)
    assert all(len(gc) == 1 for gc in mapped)

    recovered = list(strat.backward_map(tiling, mapped, children))
    assert len(recovered) == 1
    assert recovered[0] == parent_gc


def test_factor_strategy_maps_infer_children_when_omitted():
    tiling = _disjoint_two_cell_tiling()
    strat = next(FactorFactory()(tiling))
    parent_gc = GriddedChord(
        Chord((0, 0, 1, 1)),
        ((0, 0), (0, 0), (1, 1), (1, 1)),
    )
    mapped = strat.forward_map(tiling, parent_gc)
    recovered = list(strat.backward_map(tiling, mapped))
    assert recovered[0] == parent_gc


def test_factor_strategy_extra_parameters_without_assumptions():
    tiling = _disjoint_two_cell_tiling()
    strat = next(FactorFactory()(tiling))
    children = strat.decomposition_function(tiling)
    extras = strat.extra_parameters(tiling, children)
    assert extras == ({}, {})


def test_factor_strategy_extra_parameters_with_assumptions():
    ass = TrackingAssumption(
        [
            GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
            GriddedChord(Chord((0, 0)), ((1, 1), (1, 1))),
        ]
    )
    tiling = Tiling(
        obstructions=(
            crossed_chord((0, 0)),
            crossed_chord((1, 1)),
        ),
        assumptions=(ass,),
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
    )
    strat = next(FactorFactory()(tiling))
    children = strat.decomposition_function(tiling)
    extras = strat.extra_parameters(tiling, children)
    assert len(extras) == 2
    # Each child inherits the portion of the assumption supported on its cells.
    assert extras == ({"k_0": "k_0"}, {"k_0": "k_0"})


def test_factor_factory_respects_workable_and_ignore_parent_flags():
    tiling = _disjoint_two_cell_tiling()
    strat = next(
        FactorFactory(unions=False, ignore_parent=False, workable=False)(tiling)
    )
    assert strat.ignore_parent is False
    assert strat.workable is False


if __name__ == "__main__":
    test_factor_decomposes_into_multiple_factors()
    test_factor_not_factorable_when_single_component()
    test_factor_components_are_cached()
    test_factor_unites_same_row_or_column()
    test_factor_unites_cells_sharing_an_obstruction()
    test_factor_unites_cells_sharing_a_requirement()
    test_factor_unites_linkage_cells()
    test_factor_unites_cells_inside_tracked_assumption_chord()
    test_factor_factors_obs_and_reqs_split_correctly()
    test_factor_empty_tiling_obs_and_reqs()
    test_factor_reducible_factorisations_groups_irreducible_factors()
    test_factor_same_row_or_col_helper()
    test_factor_factory_yields_a_strategy()
    test_factor_factory_smoke_on_small_tiling()
    test_factor_factory_yields_nothing_when_not_factorable()
    test_factor_factory_emits_irreducible_partition()
    test_factor_factory_unions_include_grouped_partitions()
    test_factor_strategy_formal_step_and_repr()
    test_factor_strategy_json_roundtrip()
    test_factor_factory_json_roundtrip()
    test_factor_strategy_forward_backward_map_roundtrip()
    test_factor_strategy_maps_infer_children_when_omitted()
    test_factor_strategy_extra_parameters_without_assumptions()
    test_factor_strategy_extra_parameters_with_assumptions()
    test_factor_factory_respects_workable_and_ignore_parent_flags()
    print("All factor tests passed.")
