import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from tiling import Tiling
from chords import Chord, GriddedChord
from algorithms.row_col_sep import RowColSeparation
from strategies.row_col_separation import RowColumnSeparationStrategy

non_crossing = Tiling(obstructions=(GriddedChord(Chord((0,)), ((0, 1),)),
                                    GriddedChord(Chord((0,)), ((1, 0),)),
                                    GriddedChord(Chord((0,)), ((2, 1),)),
                                    GriddedChord(Chord((0,)), ((3, 0),)),
                                    GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
                                    GriddedChord(Chord((0, 0)), ((1, 1), (3, 1))),
                                    GriddedChord(Chord((0, 0)), ((2, 0), (2, 0))),
                                    GriddedChord(Chord((0, 1)), ((0, 0), (0, 0))),
                                    GriddedChord(Chord((0, 1)), ((0, 0), (2, 0))),
                                    GriddedChord(Chord((0, 1)), ((2, 0), (2, 0))),
                                    GriddedChord(Chord((0, 1, 0, 1)), ((1, 1), (1, 1), (1, 1), (1, 1))),
                                    GriddedChord(Chord((0, 1, 0, 1)), ((3, 1), (3, 1), (3, 1), (3, 1))),
                                    GriddedChord(Chord((1, 0)), ((0, 0), (0, 0))),
                                    GriddedChord(Chord((1, 0)), ((0, 0), (2, 0))),
                                    GriddedChord(Chord((1, 0)), ((1, 1), (3, 1))),
                                    GriddedChord(Chord((1, 0)), ((2, 0), (2, 0)))),
                      requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),),), 
                      simplify=False)

nc_separated = Tiling(obstructions=(GriddedChord(Chord((0,)), ((0, 1),)),
                                    GriddedChord(Chord((0,)), ((0, 2),)),
                                    GriddedChord(Chord((0,)), ((1, 0),)),
                                    GriddedChord(Chord((0,)), ((1, 1),)),
                                    GriddedChord(Chord((0,)), ((2, 1),)),
                                    GriddedChord(Chord((0,)), ((2, 2),)),
                                    GriddedChord(Chord((0,)), ((3, 0),)),
                                    GriddedChord(Chord((0,)), ((3, 2),)),
                                    GriddedChord(Chord((0, 1, 0, 1)), ((1, 2), (1, 2), (1, 2), (1, 2))),
                                    GriddedChord(Chord((0, 1, 0, 1)), ((3, 1), (3, 1), (3, 1), (3, 1)))),
                      requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (2, 0))),),))

sep = RowColSeparation(non_crossing)

assert sep.separated_tiling() == nc_separated


row_col_sep_strat = RowColumnSeparationStrategy()

assert row_col_sep_strat.decomposition_function(non_crossing) == (nc_separated,)



