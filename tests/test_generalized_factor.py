import _direct_run_bootstrap

from chord_diagrams.chords import Chord, GriddedChord
from chord_diagrams.tiling import Tiling

from chord_diagrams.algorithms.generalized_factor import GeneralizedFactor
from chord_diagrams.strategies.generalized_factor import (
    GeneralizedFactorFactory,
    GeneralizedFactorStrategy,
)


def _atom_cell_constraints(cell):
    """
    Obstructions required for building an empty cell tiling.
    (Forbid all size-2 chord patterns in the cell)
    """
    return (
        # forbid the 3 size-2 chord patterns in this cell
        GriddedChord(Chord((0, 0, 1, 1)), (cell,) * 4),
        GriddedChord(Chord((0, 1, 1, 0)), (cell,) * 4),
        GriddedChord(Chord((0, 1, 0, 1)), (cell,) * 4),
    )


def _noncrossing_cell_obstruction(cell):
    """Noncrossing chord pattern in the cell"""
    return GriddedChord(Chord((0, 1, 0, 1)), (cell,) * 4)


def test_generalized_factor_finds_shared_atom_cover():
    # Arrange: three active cells on distinct rows/cols.
    u = (0, 0)  # shared overlap U (atom)
    a = (1, 2)
    b = (3, 1)

    tiling = Tiling(
        obstructions=(
            *_atom_cell_constraints(u),
            _noncrossing_cell_obstruction(a),
            _noncrossing_cell_obstruction(b),
        ),
        requirements=((GriddedChord(Chord((0, 0)), (u, u)),),),
        simplify=False,
        remove_empty_rows_and_cols=False,
        expand=False,
    )

    # Sanity: the overlap subtiling should be an atom.
    u_tiling = tiling.sub_tiling((u,))
    assert u_tiling.is_atom()
    n_u = u_tiling.minimum_size_of_object()

    print(tiling)
    #tiling.pretty_print_latex("tiling.tex")
    print("tiling printed")

    # Act
    gf = GeneralizedFactor(tiling)
    results = list(gf.generalized_factorizations())
    print("results: ", results)

    # Assert
    assert len(results) == 1
    parts, shifts, shift_rel = results[0]

    # We expect two parts, each sharing the overlap cell u.
    assert len(parts) == 2
    assert set(parts[0]) == {a, u} or set(parts[1]) == {a, u}
    assert set(parts[0]) == {b, u} or set(parts[1]) == {b, u}

    # Shifts should correspond to the overlap size (up to ordering).
    assert shifts.count(0) == 1
    assert shifts.count(n_u) == 1

    # shift_reliance is allowed to be empty here; just ensure it's a set of indices.
    assert isinstance(shift_rel, set)
    assert all(isinstance(i, int) for i in shift_rel)


def test_generalized_factor_factory_emits_strategy_with_shifts():
    u = (0, 0)
    a = (1, 2)
    b = (3, 1)
    tiling = Tiling(
        obstructions=(
            *_atom_cell_constraints(u),
            _noncrossing_cell_obstruction(a),
            _noncrossing_cell_obstruction(b),
        ),
        requirements=((GriddedChord(Chord((0, 0)), (u, u)),),),
        simplify=False,
        remove_empty_rows_and_cols=False,
        expand=False,
    )

    factory = GeneralizedFactorFactory()
    strats = list(factory(tiling))
    assert strats, "Expected GeneralizedFactorFactory to yield a strategy"
    strat = strats[0]
    assert isinstance(strat, GeneralizedFactorStrategy)

    children = strat.decomposition_function(tiling)
    assert len(children) == len(strat.partition)
    assert len(strat.shifts) == len(strat.partition)


def test_generalized_factor_searcher_small_example():
    """
    Use CombinatorialSpecificationSearcher (as in example.py) on a small tiling
    where generalized factorization exists, and ensure a specification can be found.
    """
    from comb_spec_searcher import AtomStrategy, CombinatorialSpecificationSearcher, StrategyPack

    u = (0, 0)  # overlap U (atom)
    a = (1, 2)
    b = (3, 1)
    tiling = Tiling(
        obstructions=(
            *_atom_cell_constraints(u),
            *_atom_cell_constraints(a),
            *_atom_cell_constraints(b),
            _noncrossing_cell_obstruction(a),
            _noncrossing_cell_obstruction(b),
        ),
        # Make each single-cell subtiling an atom by requiring a 00 in each.
        requirements=(
            (GriddedChord(Chord((0, 0)), (u, u)),),
            (GriddedChord(Chord((0, 0)), (a, a)),),
            (GriddedChord(Chord((0, 0)), (b, b)),),
        ),
        simplify=False,
        remove_empty_rows_and_cols=False,
        expand=False,
    )

    # Strategy pack: keep it minimal. GeneralizedFactorFactory prefers generalized
    # factorization (nonzero shifts) when it exists, then falls back to disjoint
    # factorization on the children, which will be single-cell atoms here.
    pack = StrategyPack(
        initial_strats=[GeneralizedFactorFactory()],
        inferral_strats=[],
        expansion_strats=[],
        ver_strats=[AtomStrategy()],
        name=("Generalized factor searcher smoke test"),
    )

    # Ensure generalized factorization is available at the root (nonzero shifts).
    root_strats = list(GeneralizedFactorFactory()(tiling))
    assert root_strats
    assert any(root_strats[0].shifts), "Expected generalized factorization (nonzero shifts) at root"

    searcher = CombinatorialSpecificationSearcher(tiling, pack)
    spec = searcher.auto_search()
    assert spec is not None
    assert spec.get_genf() is not None


test_generalized_factor_finds_shared_atom_cover()
test_generalized_factor_factory_emits_strategy_with_shifts()
test_generalized_factor_searcher_small_example()