import _direct_run_bootstrap 

from src.common.chords import Chord, GriddedChord
from src.common.obstructions import (
    FiniteObstruction,
    InfiniteObstruction,
    InfinitePatternType,
    Obstruction,
)


def _single_cell(chord: Chord, cell=(0, 0)) -> GriddedChord:
    return GriddedChord.single_cell(chord, cell)


# --- Chord pattern constructors ---

assert Chord.top_cycle(3) == Chord.bottom_cycle(3) == Chord((0, 1, 2, 0, 1, 2))
assert Chord.top_cycle(4) == Chord((0, 1, 2, 0, 3, 2, 1, 3))
assert Chord.bottom_cycle(4) == Chord((0, 1, 2, 3, 1, 0, 3, 2))
assert Chord.partial_top_cycle(3) == Chord((0, 1, 2, 1, 0, 2))
assert Chord.partial_bottom_cycle(3) == Chord((0, 1, 2, 0, 2, 1))
assert Chord.nonnesting_path(2) == Chord((0, 1, 0, 1))
assert Chord.nonnesting_path(3) == Chord((0, 1, 0, 2, 1, 2))
assert Chord.top_cycle(5).remove_chord(0) == Chord.partial_top_cycle(4)
assert Chord.bottom_cycle(5).remove_chord(0) == Chord.partial_bottom_cycle(4)


# --- FiniteObstruction ---

gc_cross = _single_cell(Chord((0, 1, 0, 1)))
gc_nest = _single_cell(Chord((0, 0, 1, 1)))
gc_triangle = _single_cell(Chord.top_cycle(3))

fin = FiniteObstruction(gc_cross)
assert fin.occurs_in(gc_triangle)
assert not fin.occurs_in(gc_nest)
assert fin.is_avoided_by(gc_nest)
assert fin.patt == gc_cross.patt
assert fin.gc == gc_cross
assert fin.cells == gc_cross.get_active_cells()

fin2 = FiniteObstruction.from_gridded_chord(gc_cross)
assert fin == fin2
assert hash(fin) == hash(fin2)
assert Obstruction.from_dict(fin.to_jsonable()) == fin


# --- InfiniteObstruction: cycles in one cell ---

top_obs = InfiniteObstruction(InfinitePatternType.TOP_CYCLES, ((0, 0),))
bottom_obs = InfiniteObstruction(InfinitePatternType.BOTTOM_CYCLES, ((0, 0),))
cycles_obs = InfiniteObstruction(InfinitePatternType.CYCLES, ((0, 0),))

assert top_obs.occurs_in(gc_triangle)
assert bottom_obs.occurs_in(gc_triangle)
assert cycles_obs.occurs_in(gc_triangle)
assert not top_obs.occurs_in(gc_cross)
assert not top_obs.occurs_in(gc_nest)

gc_top4 = _single_cell(Chord.top_cycle(4))
gc_bottom4 = _single_cell(Chord.bottom_cycle(4))
assert top_obs.occurs_in(gc_top4)
assert not bottom_obs.occurs_in(gc_top4)
assert bottom_obs.occurs_in(gc_bottom4)
assert not top_obs.occurs_in(gc_bottom4)
assert cycles_obs.occurs_in(gc_top4)
assert cycles_obs.occurs_in(gc_bottom4)


# --- Partial cycles ---

partial_top_obs = InfiniteObstruction(
    InfinitePatternType.PARTIAL_TOP_CYCLES, ((0, 0),)
)
partial_bottom_obs = InfiniteObstruction(
    InfinitePatternType.PARTIAL_BOTTOM_CYCLES, ((0, 0),)
)
partial_cycles_obs = InfiniteObstruction(
    InfinitePatternType.PARTIAL_CYCLES, ((0, 0),)
)

gc_pt3 = _single_cell(Chord.partial_top_cycle(3))
gc_pb3 = _single_cell(Chord.partial_bottom_cycle(3))
assert partial_top_obs.occurs_in(gc_pt3)
assert not partial_top_obs.occurs_in(gc_pb3)
assert partial_bottom_obs.occurs_in(gc_pb3)
assert not partial_bottom_obs.occurs_in(gc_pt3)
assert partial_cycles_obs.occurs_in(gc_pt3)
assert partial_cycles_obs.occurs_in(gc_pb3)


# --- Nonnesting paths ---

path_obs = InfiniteObstruction(InfinitePatternType.NONNESTING_PATHS, ((0, 0),))
gc_path2 = _single_cell(Chord.nonnesting_path(2))
gc_path3 = _single_cell(Chord.nonnesting_path(3))
assert path_obs.occurs_in(gc_path2)
assert path_obs.occurs_in(gc_path3)
assert path_obs.occurs_in(gc_triangle)  # triangle contains a path of length 2
assert not path_obs.is_avoided_by(gc_path2)
assert path_obs.is_avoided_by(gc_nest)


# --- Cell-boundary constraints ---

# Triangle entirely in (1, 0) should not match an obstruction that requires
# first and last points in (0, 0).
obs_two_cells = InfiniteObstruction(
    InfinitePatternType.TOP_CYCLES, ((0, 0), (1, 0))
)
gc_triangle_right = _single_cell(Chord.top_cycle(3), (1, 0))
assert not obs_two_cells.occurs_in(gc_triangle_right)
assert not obs_two_cells.occurs_in(gc_triangle)  # all points in (0,0), last != (1,0)

# A top cycle gridded with first point in (0,0) and last in (1,0).
# Pattern T3 = (0,1,2,0,1,2): assign cells so first is e1 and last is ek.
positions = ((0, 0), (0, 0), (0, 0), (1, 0), (1, 0), (1, 0))
gc_split = GriddedChord(Chord.top_cycle(3), positions)
assert obs_two_cells.occurs_in(gc_split)

# Wrong start cell: first point not in e1.
bad_positions = ((1, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0))
gc_bad_start = GriddedChord(Chord.top_cycle(3), bad_positions)
assert not obs_two_cells.occurs_in(gc_bad_start)

# Partial top: second point must also be in e1.
partial_obs_cells = InfiniteObstruction(
    InfinitePatternType.PARTIAL_TOP_CYCLES, ((0, 0), (1, 0))
)
# partial top 3 = (0,1,2,1,0,2); require pos[0]=pos[1]=(0,0), pos[-1]=(1,0)
pt_ok = GriddedChord(
    Chord.partial_top_cycle(3),
    ((0, 0), (0, 0), (0, 0), (0, 0), (1, 0), (1, 0)),
)
pt_bad_second = GriddedChord(
    Chord.partial_top_cycle(3),
    ((0, 0), (1, 0), (1, 0), (1, 0), (1, 0), (1, 0)),
)
assert partial_obs_cells.occurs_in(pt_ok)
assert not partial_obs_cells.occurs_in(pt_bad_second)


# --- JSON / equality / sorting ---

inf = InfiniteObstruction(InfinitePatternType.TOP_CYCLES, ((1, 0), (0, 0)))
# cells are canonicalized to lex order
assert inf.pos == ((0, 0), (1, 0))
assert Obstruction.from_dict(inf.to_jsonable()) == inf

assert fin < inf  # finite sorts before infinite
assert sorted([inf, fin]) == [fin, inf]

# string pattern_type accepted
inf2 = InfiniteObstruction("top_cycles", ((0, 0),))
assert inf2 == top_obs

# --- Tiling integration ---

from src.common.tiling import Tiling

t_finite = Tiling((gc_cross,), simplify=False, expand=False, remove_empty_rows_and_cols=False)
assert len(t_finite.obstructions) == 1
assert isinstance(t_finite.obstructions[0], FiniteObstruction)
assert t_finite.obstructions[0].gc == gc_cross

t_inf = Tiling(
    (InfiniteObstruction(InfinitePatternType.TOP_CYCLES, ((0, 0),)),),
    simplify=False,
    expand=False,
    remove_empty_rows_and_cols=False,
)
assert isinstance(t_inf.obstructions[0], InfiniteObstruction)
assert not t_inf.contains(gc_triangle)
assert t_inf.contains(gc_nest)

assert Tiling.from_dict(t_inf.to_jsonable()) == t_inf

print("test_obstructions: all checks passed")

