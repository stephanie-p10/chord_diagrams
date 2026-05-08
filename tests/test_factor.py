import _direct_run_bootstrap 

from chord_diagrams.algorithms.factor import Factor
from chord_diagrams.chords import Chord, GriddedChord
from chord_diagrams.strategies.factor import FactorFactory, FactorStrategy
from chord_diagrams.tiling import Tiling


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


def test_factor_decomposes_into_multiple_factors():
    tiling = _make_non_crossing()
    algo = Factor(tiling)
    assert algo.factorable()

    factors = algo.factors()
    assert isinstance(factors, tuple)
    assert len(factors) >= 2
    assert all(isinstance(t, Tiling) for t in factors)


def test_factor_factory_yields_a_strategy():
    tiling = _make_non_crossing()
    factory = FactorFactory(True)
    strats = list(factory(tiling))
    assert strats
    assert isinstance(strats[0], FactorStrategy)
