import _direct_run_bootstrap

from src.algorithms.chord_placement import (
    ChordPlacement,
    RequirementPlacement,
    stretch_cells_around_placement,
    stretch_infinite_obstruction,
)
from src.common import DIR_SOUTH
from src.common.chords import Chord, GriddedChord
from src.common.obstructions import InfiniteObstruction, InfinitePatternType
from src.common.tiling import Tiling


# --- cell stretch helper (mirrors linkage mu_e) ---

assert stretch_cells_around_placement(((0, 0),), (0, 0), True, True) == (
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 0),
    (1, 1),
    (1, 2),
    (2, 0),
    (2, 1),
    (2, 2),
)
assert stretch_cells_around_placement(((1, 0),), (0, 0), True, True) == ((3, 0), (3, 1), (3, 2))
assert stretch_cells_around_placement(((0, 0),), (1, 0), True, True) == ((0, 0), (0, 1), (0, 2))
assert stretch_cells_around_placement(((0, 0),), (0, 0), False, True) == (
    (0, 0),
    (1, 0),
    (2, 0),
)

inf = InfiniteObstruction(InfinitePatternType.TOP_CYCLES, ((0, 0),))
stretched = stretch_infinite_obstruction(inf, (0, 0), True, True)
assert len(stretched) == 9
assert all(ob.pattern_type == InfinitePatternType.TOP_CYCLES for ob in stretched)
assert {ob.pos for ob in stretched} == {
    ((0, 0),),
    ((0, 1),),
    ((0, 2),),
    ((1, 0),),
    ((1, 1),),
    ((1, 2),),
    ((2, 0),),
    ((2, 1),),
    ((2, 2),),
}

# Multi-cell family: unique shift of the non-placed cell, product with images of placed.
inf2 = InfiniteObstruction(InfinitePatternType.BOTTOM_CYCLES, ((0, 0), (1, 0)))
stretched2 = stretch_infinite_obstruction(inf2, (0, 0), True, True)
assert all(len(ob.pos) == 2 for ob in stretched2)
assert all(ob.pos[-1][0] == 3 for ob in stretched2)  # (1,0) shifts to x=3


# --- multiplex helpers keep infinite obs ---

tiling = Tiling(
    obstructions=(
        InfiniteObstruction(InfinitePatternType.TOP_CYCLES, ((0, 0),)),
        GriddedChord(Chord((0, 1, 0, 1)), ((0, 0),) * 4),
    ),
    requirements=((GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),),),
    simplify=False,
    expand=False,
    remove_empty_rows_and_cols=False,
    derive_empty=False,
)

req_pl = RequirementPlacement(tiling)
stretched_obs = list(req_pl.stretched_obs((0, 0)))
inf_stretched = [ob for ob in stretched_obs if isinstance(ob, InfiniteObstruction)]
assert len(inf_stretched) == 9

chord_pl = ChordPlacement(tiling)
multiplexed = chord_pl.get_multiplexes_of_gcs(tiling.obstructions, (0, 0))
inf_mux = [ob for ob in multiplexed if isinstance(ob, InfiniteObstruction)]
assert len(inf_mux) == 9
assert set(inf_mux) == set(inf_stretched)


# --- placement must retain remapped single-cell infinite obstructions ---

placed = chord_pl.place_chord(
    GriddedChord(Chord((0, 0)), ((0, 0), (0, 0))),
    0,
    DIR_SOUTH,
)
inf_after = [ob for ob in placed.obstructions if isinstance(ob, InfiniteObstruction)]
assert len(inf_after) >= 1
assert all(ob.pattern_type == InfinitePatternType.TOP_CYCLES for ob in inf_after)

# Some surviving single-cell family still forbids a top cycle in that cell.
single_cell_inf = [ob for ob in inf_after if len(ob.pos) == 1]
assert single_cell_inf
found = False
for survivor in single_cell_inf:
    cycle_cell = survivor.pos[0]
    host = GriddedChord.single_cell(Chord.top_cycle(3), cycle_cell)
    if survivor.occurs_in(host):
        assert not placed.contains(host)
        found = True
        break
assert found, "expected a surviving single-cell top-cycle obstruction"

print("test_infinite_placement: ok")
