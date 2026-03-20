from tiling import Tiling, Chord, GriddedChord
from collections import Counter
from strategies.factor import FactorFactory, FactorStrategy, Factor
import time

from strategies.obstruction_inferral import SubobstructionInferralFactory


# first pattern in theorem 3.1.2 class
c1 = Chord((0, 1, 2, 0, 1, 2))
# second pattern in theorem 3.1.2 class
c2 = Chord((0, 1, 2, 0, 2, 1))
single_chord = Chord((0, 0))

'''my_t3_prime = Tiling(obstructions=(GriddedChord.single_chord(((0, 0), (0, 0))),
                                GriddedChord.single_chord(((2, 0), (2, 0))), 
                                GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0))),

                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),

                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),

                                GriddedChord(c1, ((0, 0), (1, 1), (1, 1), (2, 0), (3, 1), (3, 1))),

                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (3, 1), (3, 1), (3, 1), (3, 1))),

                                GriddedChord(c2, ((0, 0), (1, 1), (1, 1), (2, 0), (3, 1), (3, 1))),

                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (3, 1), (3, 1), (3, 1), (3, 1))),),
                  requirements=((GriddedChord.single_chord(((0, 0), (2, 0))),),))

lukas_t3_prime = Tiling(obstructions=(GriddedChord.single_chord(((0, 0), (0, 0))),
                                GriddedChord.single_chord(((2, 0), (2, 0))), 
                                GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (2, 0), (2, 0))),
                                GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (2, 0), (2, 0))),

                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                GriddedChord(c1, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                GriddedChord(c1, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c1, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),

                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))),
                                GriddedChord(c2, ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (3, 1))),
                                GriddedChord(c2, ((1, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                GriddedChord(c2, ((3, 1), (3, 1), (3, 1), (3, 1), (3, 1), (3, 1))),
                                
                                GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (3, 1), (3, 1))),
                                GriddedChord(Chord((0, 1, 1, 0)), ((1, 1), (1, 1), (3, 1), (3, 1))),),
                  requirements=((GriddedChord.single_chord(((0, 0), (2, 0))),),))

t6 = Tiling(obstructions=(GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                          GriddedChord(Chord((0, 0)), ((1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (0, 0), (1, 0), (1, 0))),
                          GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (0, 0), (1, 0), (1, 0))),),
            requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (1, 0))),),))'''

'''non_crossing_ex1 = Tiling((GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (1, 1), (1, 1))),
                       GriddedChord(Chord((0, 1, 0, 1)), ((3, 1), (3, 1), (3, 1), (3, 1))), 
                       GriddedChord(Chord((0, 0)), ((1, 1), (3, 1))), 
                       GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))), 

                       GriddedChord(Chord((0, 1)), ((1, 1), (3, 1))), 
                       GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))), 
                       GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((1, 0)), ((2, 0), (2, 0))), ),
                       ((GriddedChord(Chord((0, 0)), ((0, 0), (3, 0))),),))'''

non_crossing_ex2 = Tiling((GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 1), (2, 0), (3, 1))),

                       GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))), 
                       GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))), 
                       GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))), 
                       GriddedChord(Chord((1, 0)), ((2, 0), (2, 0))), ),
                       ((GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),),))


factory = SubobstructionInferralFactory()

for strat in factory(non_crossing_ex2):
    add_subobs_tiling = strat.decomposition_function(non_crossing_ex2)[0]
    add_subobs_tiling._simplify()
    print(add_subobs_tiling)



