import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from chords import Chord, GriddedChord
from tiling import Tiling
from algorithms.factor import Factor
from strategies.factor import FactorFactory, FactorStrategy

non_crossing = Tiling(obstructions=sorted((GriddedChord(Chord((0,)), ((0, 1),)),
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
                                    #GriddedChord(Chord((1, 0)), ((1, 2), (3, 1))), # it should have gotten rid of this...
                                    GriddedChord(Chord((0, 1, 0, 1)), ((1, 2), (1, 2), (1, 2), (1, 2))),
                                    GriddedChord(Chord((0, 1, 0, 1)), ((3, 1), (3, 1), (3, 1), (3, 1)))),),
                      requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),),))


#for ob in non_crossing.obstructions:
    #print(ob)
#print(non_crossing.active_cells)

factor_nc_algo = Factor(non_crossing)
factors = factor_nc_algo.factors()

nc_factors_list = [Tiling(obstructions=(GriddedChord(Chord((0,1,0,1)), ((0,0), (0,0), (0,0), (0,0))),)),
                    Tiling(obstructions=(GriddedChord(Chord((0,1,0,1)), ((0,0), (0,0), (0,0), (0,0))),)),
                    Tiling(obstructions=(GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                                         GriddedChord(Chord((0, 0)), ((1, 0), (1, 0))),
                                         GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))),
                                         GriddedChord(Chord((0, 1)), ((0, 0), (1, 0))),
                                         GriddedChord(Chord((0, 1)), ((1, 0), (1, 0))),
                                         GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))),
                                         GriddedChord(Chord((1, 0)), ((0, 0), (1, 0))),
                                         GriddedChord(Chord((1, 0)), ((1, 0), (1, 0)))),
                                         requirements=((GriddedChord(Chord((0,0)), ((0, 0), (1, 0))),),))]

for factor_list in [factor for factor in factor_nc_algo.factors()]:
    print(factor_list)
assert [factor for factor in factor_nc_algo.factors()] == nc_factors_list


factory = FactorFactory(True)
for strat in factory(non_crossing):
    print("next strat:")
    for part in strat.decomposition_function(non_crossing):
        print(part)
    print()
