import _direct_run_bootstrap  

from chord_diagrams.algorithms.obstruction_inferral import (
    AllObstructionInferral,
    EmptyCellInferral,
    SubobstructionInferral,
)
from chord_diagrams.chords import Chord, GriddedChord
from chord_diagrams.strategies.obstruction_inferral import ObstructionInferralFactory
from chord_diagrams.tiling import Tiling


def crossed_chord(pos):
    return GriddedChord(Chord((0, 1, 0, 1)), (pos,) * 4)


def single_chord(pos_left, pos_right):
    return GriddedChord(Chord((0, 0)), (pos_left, pos_right))


def _make_non_crossing() -> Tiling:
    return Tiling(
        obstructions=(
            single_chord((1, 1), (3, 1)),
            crossed_chord((1, 1)),
            crossed_chord((3, 1)),
            single_chord((0, 0), (0, 0)),
            single_chord((2, 0), (2, 0)),
        ),
        requirements=([single_chord((0, 0), (2, 0))],),
        simplify=False,
        expand=False,
    )

def test_obstruction_inferral_algorithms_smoke():
    tiling = _make_non_crossing()
    assert isinstance(AllObstructionInferral(tiling, 2).new_obs(), list)
    assert isinstance(SubobstructionInferral(tiling).new_obs(), list)
    assert isinstance(EmptyCellInferral(tiling).new_obs(), list)


def test_obstruction_inferral_factory_smoke():
    tiling = _make_non_crossing()
    factory = ObstructionInferralFactory()
    # May yield or not depending on tiling; just ensure no crash.
    list(factory(tiling))

