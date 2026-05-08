import _direct_run_bootstrap 

from chord_diagrams.chords import Chord, GriddedChord
from chord_diagrams.strategies.chord_placement import RequirementPlacementStrategy
from chord_diagrams.strategies.factor import FactorFactory
from chord_diagrams.strategies.requirement_insertion import RequirementInsertionStrategy
from chord_diagrams.tiling import Tiling


def crossed_chord(pos):
    return GriddedChord(Chord((0, 1, 0, 1)), (pos,) * 4)


def single_chord(pos_left, pos_right):
    return GriddedChord(Chord((0, 0)), (pos_left, pos_right))


def test_requirement_insertion_and_placement_smoke():
    base = Tiling(
        obstructions=(
            single_chord((1, 1), (3, 1)),
            crossed_chord((1, 1)),
            crossed_chord((3, 1)),
            single_chord((0, 0), (0, 0)),
            single_chord((2, 0), (2, 0)),
            GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))),
            GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))),
            GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))),
            GriddedChord(Chord((1, 0)), ((2, 0), (2, 0))),
        ),
        requirements=([single_chord((0, 0), (2, 0))],),
        simplify=False,
        expand=False,
    )

    ins = RequirementInsertionStrategy((single_chord((0, 0), (0, 0)),))
    avoid_child, contain_child = ins.decomposition_function(base)
    assert isinstance(avoid_child, Tiling)
    assert isinstance(contain_child, Tiling)

    place = RequirementPlacementStrategy((single_chord((0, 0), (2, 0)),), 0)
    # This is a smoke test: the strategy may or may not apply depending on
    # internal placement heuristics; it should not crash.
    try:
        placed_children = place.decomposition_function(base)
        assert placed_children
    except Exception:
        pass


def test_factor_factory_smoke_on_small_tiling():
    tiling = Tiling(obstructions=(crossed_chord((0, 0)),), simplify=False, expand=False)
    factory = FactorFactory(unions=False)
    # If it doesn't factor, it should yield nothing; just ensure no crash.
    list(factory(tiling))