"""Obstruction types for tilings of gridded chord diagrams.

Finite obstructions wrap a single :class:`~src.common.chords.GriddedChord`.
Infinite obstructions compactly represent Nabergall's infinite families
(top/bottom cycles, partials, and nonnesting paths) constrained to a list of
cells.

``Tiling`` stores a tuple of :class:`Obstruction`. Bare ``GriddedChord`` values
passed to ``Tiling`` are normalized to :class:`FiniteObstruction`.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from itertools import combinations
from typing import FrozenSet, Iterable, Iterator, Sequence, Tuple, Union

from src.common.chords import Chord, GriddedChord

Cell = Tuple[int, int]

ObstructionLike = Union["Obstruction", GriddedChord]


def normalize_obstruction(ob: ObstructionLike) -> "Obstruction":
    """Coerce a ``GriddedChord`` or ``Obstruction`` to an ``Obstruction``."""
    if isinstance(ob, Obstruction):
        return ob
    if isinstance(ob, GriddedChord):
        return FiniteObstruction(ob)
    raise TypeError(f"Expected Obstruction or GriddedChord, got {type(ob)}")


def normalize_obstructions(obs: Iterable[ObstructionLike]) -> Tuple["Obstruction", ...]:
    return tuple(normalize_obstruction(ob) for ob in obs)


class InfinitePatternType(Enum):
    # Possible types of infinite patterns

    CYCLES = "cycles"  # G_{{T^{>=3}, B^{>=3}}}
    TOP_CYCLES = "top_cycles"  # G_{T^{>=3}}
    BOTTOM_CYCLES = "bottom_cycles"  # G_{B^{>=3}}
    PARTIAL_TOP_CYCLES = "partial_top_cycles"  # G_{T'^{>=3}}
    PARTIAL_BOTTOM_CYCLES = "partial_bottom_cycles"  # G_{B'^{>=3}}
    PARTIAL_CYCLES = "partial_cycles"  # G_{{T'^{>=3}, B'^{>=3}}}
    NONNESTING_PATHS = "nonnesting_paths"  # G_{P^{>1}}


class Obstruction(ABC):
    # Obstructions base class

    @property
    @abstractmethod
    def cells(self) -> FrozenSet[Cell]:
        ...

    @property
    def pos(self) -> Tuple[Cell, ...]:
        """Ordered cell data used by callers that expect a ``.pos`` attribute."""
        return self._pos_tuple()

    @abstractmethod
    def _pos_tuple(self) -> Tuple[Cell, ...]:     
        ...

    @abstractmethod
    def occurs_in(self, gc: GriddedChord) -> bool:
        ...

    def is_avoided_by(self, gc: GriddedChord) -> bool:
        return not self.occurs_in(gc)

    def is_point(self) -> bool:
        return False

    def is_empty(self) -> bool:
        return False

    def is_localized(self) -> bool:
        return False

    def is_single_cell(self) -> bool:
        return False

    def is_single_chord(self) -> bool:
        return False

    @abstractmethod
    def occurrences_in(self, other: GriddedChord) -> Iterator:
        """Compatibility with ``GriddedChord.contains`` / ``avoids``."""

    @abstractmethod
    def _sort_key(self) -> tuple:
        ...

    @abstractmethod
    def to_jsonable(self) -> dict:
        ...

    @classmethod
    def from_dict(cls, d: dict) -> "Obstruction":
        kind = d.get("kind")
        if kind == "finite":
            return FiniteObstruction.from_dict(d)
        if kind == "infinite":
            return InfiniteObstruction.from_dict(d)
        # Backward compatibility: raw gridded-chord JSON without a kind field.
        if "patt" in d and "pos" in d:
            return FiniteObstruction(GriddedChord.from_dict(d))
        raise ValueError(f"Unknown obstruction kind: {kind!r}")

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, Obstruction):
            return NotImplemented
        return self._sort_key() < other._sort_key()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Obstruction):
            return NotImplemented
        return self._sort_key() == other._sort_key()

    def __hash__(self) -> int:
        return hash(self._sort_key())


class FiniteObstruction(Obstruction):
    """A finite obstruction: a single gridded chord pattern to avoid."""

    def __init__(self, gc: GriddedChord) -> None:
        if not isinstance(gc, GriddedChord):
            raise TypeError(f"FiniteObstruction requires a GriddedChord, got {type(gc)}")
        self._gc = gc

    @classmethod
    def from_gridded_chord(cls, gc: GriddedChord) -> "FiniteObstruction":
        return cls(gc)

    @property
    def gc(self) -> GriddedChord:
        return self._gc

    @property
    def patt(self) -> Tuple[int, ...]:
        return self._gc.patt

    @property
    def cells(self) -> FrozenSet[Cell]:
        return self._gc.get_active_cells()

    def _pos_tuple(self) -> Tuple[Cell, ...]:
        return self._gc.pos

    def occurs_in(self, gc: GriddedChord) -> bool:
        return self._gc.occurs_in(gc)

    def occurrences_in(self, other: GriddedChord) -> Iterator:
        return self._gc.occurrences_in(other)

    def is_point(self) -> bool:
        return self._gc.is_point()

    def is_empty(self) -> bool:
        return self._gc.is_empty()

    def is_localized(self) -> bool:
        return self._gc.is_localized()

    def is_single_cell(self) -> bool:
        return self._gc.is_single_cell()

    def is_single_chord(self) -> bool:
        return self._gc.is_single_chord()

    @property
    def _cells(self) -> FrozenSet[Cell]:
        return self._gc._cells

    def __len__(self) -> int:
        return len(self._gc)

    def __getattr__(self, name: str):
        # Forward remaining GriddedChord APIs (e.g. all_subchords, contradictory).
        return getattr(self._gc, name)

    def _sort_key(self) -> tuple:
        return (0, self._gc.patt, self._gc.pos)

    def to_jsonable(self) -> dict:
        return {"kind": "finite", "gc": self._gc.to_jsonable()}

    @classmethod
    def from_dict(cls, d: dict) -> "FiniteObstruction":
        if "gc" in d:
            return cls(GriddedChord.from_dict(d["gc"]))
        return cls(GriddedChord.from_dict(d))

    def __repr__(self) -> str:
        return f"FiniteObstruction({self._gc!r})"

    def __str__(self) -> str:
        return str(self._gc)


class InfiniteObstruction(Obstruction):
    """An infinite family of gridded patterns constrained to a cell list.

    Following Nabergall, for cells ``e1 < ... < ek`` (lex order) the family
    ``H_{e1,...,ek}`` consists of all patterns in the family whose points lie
    only in those cells, with the first point in ``e1`` and the last in ``ek``.
    Partial top (resp. bottom) families additionally require the second
    (resp. second-to-last) point to lie in ``e1`` (resp. ``ek``).
    """

    def __init__(
        self,
        pattern_type: Union[InfinitePatternType, str],
        cells: Iterable[Cell],
    ) -> None:
        if isinstance(pattern_type, str):
            pattern_type = InfinitePatternType(pattern_type)
        if not isinstance(pattern_type, InfinitePatternType):
            raise TypeError(
                f"pattern_type must be InfinitePatternType, got {type(pattern_type)}"
            )
        cell_list = tuple((int(x), int(y)) for x, y in cells)
        if not cell_list:
            raise ValueError("InfiniteObstruction requires at least one cell")
        if len(set(cell_list)) != len(cell_list):
            raise ValueError("InfiniteObstruction cells must be unique")
        ordered = tuple(sorted(cell_list))
        if ordered != cell_list:
            # Canonicalize to lexicographic order as in the paper.
            cell_list = ordered
        self._pattern_type = pattern_type
        self._cells = cell_list

    @property
    def pattern_type(self) -> InfinitePatternType:
        return self._pattern_type

    @property
    def cells(self) -> FrozenSet[Cell]:
        return frozenset(self._cells)

    def _pos_tuple(self) -> Tuple[Cell, ...]:
        return self._cells

    @property
    def _pos(self) -> Tuple[Cell, ...]:
        """Alias for callers that historically used ``ob._pos``."""
        return self._cells

    def is_single_cell(self) -> bool:
        return len(self._cells) == 1

    def occurrences_in(self, other: GriddedChord) -> Iterator:
        if self.occurs_in(other):
            yield ()

    def __len__(self) -> int:
        return 10**9

    def _sort_key(self) -> tuple:
        return (1, self._pattern_type.value, self._cells)

    def to_jsonable(self) -> dict:
        return {
            "kind": "infinite",
            "pattern_type": self._pattern_type.value,
            "cells": [list(c) for c in self._cells],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "InfiniteObstruction":
        return cls(
            InfinitePatternType(d["pattern_type"]),
            tuple(map(tuple, d["cells"])),
        )

    def __repr__(self) -> str:
        return (
            f"InfiniteObstruction({self._pattern_type!r}, {self._cells!r})"
        )

    def __str__(self) -> str:
        return f"InfiniteObstruction({self._pattern_type.value}, {self._cells})"

    def map_cells(self, cell_map) -> "InfiniteObstruction":
        """Return a copy with cells transformed by ``cell_map``."""
        return InfiniteObstruction(
            self._pattern_type,
            tuple(cell_map(cell) for cell in self._cells),
        )

    # ------------------------------------------------------------------
    # Pattern generators for each family
    # ------------------------------------------------------------------

    def _candidate_patterns(self, max_size: int) -> Iterator[Tuple[Chord, str]]:
        """Yield ``(pattern, variant)`` pairs up to ``max_size`` chords.

        ``variant`` is ``"top"``, ``"bottom"``, ``"partial_top"``,
        ``"partial_bottom"``, or ``"path"`` and controls endpoint constraints.
        """
        pt = self._pattern_type
        if pt in (
            InfinitePatternType.TOP_CYCLES,
            InfinitePatternType.BOTTOM_CYCLES,
            InfinitePatternType.CYCLES,
        ):
            for n in range(3, max_size + 1):
                if pt in (
                    InfinitePatternType.TOP_CYCLES,
                    InfinitePatternType.CYCLES,
                ):
                    yield Chord.top_cycle(n), "top"
                if pt in (
                    InfinitePatternType.BOTTOM_CYCLES,
                    InfinitePatternType.CYCLES,
                ):
                    # T_3 = B_3; skip duplicate for the union family.
                    if pt == InfinitePatternType.CYCLES and n == 3:
                        continue
                    yield Chord.bottom_cycle(n), "bottom"
        elif pt in (
            InfinitePatternType.PARTIAL_TOP_CYCLES,
            InfinitePatternType.PARTIAL_BOTTOM_CYCLES,
            InfinitePatternType.PARTIAL_CYCLES,
        ):
            for n in range(3, max_size + 1):
                if pt in (
                    InfinitePatternType.PARTIAL_TOP_CYCLES,
                    InfinitePatternType.PARTIAL_CYCLES,
                ):
                    yield Chord.partial_top_cycle(n), "partial_top"
                if pt in (
                    InfinitePatternType.PARTIAL_BOTTOM_CYCLES,
                    InfinitePatternType.PARTIAL_CYCLES,
                ):
                    yield Chord.partial_bottom_cycle(n), "partial_bottom"
        elif pt == InfinitePatternType.NONNESTING_PATHS:
            for n in range(2, max_size + 1):
                yield Chord.nonnesting_path(n), "path"
        else:
            raise ValueError(f"Unhandled pattern type: {pt}")

    @staticmethod
    def _point_occurrences(
        pattern: Chord, host: GriddedChord
    ) -> Iterator[Tuple[int, ...]]:
        """Yield increasing tuples of point indices of ``pattern`` in ``host``."""
        host_patt = host.patt
        plen = len(pattern.get_pattern())
        if plen > len(host_patt):
            return
        target = tuple(Chord.reindex(pattern.get_pattern()))
        indexed = list(enumerate(host_patt))
        for subslice in combinations(indexed, plen):
            vals = [val for _, val in subslice]
            if tuple(Chord.reindex(vals)) == target:
                yield tuple(idx for idx, _ in subslice)

    def _respects_cell_constraints(
        self, host: GriddedChord, point_indices: Sequence[int], variant: str
    ) -> bool:
        positions = [host.pos[i] for i in point_indices]
        allowed = set(self._cells)
        if any(p not in allowed for p in positions):
            return False
        if positions[0] != self._cells[0]:
            return False
        if positions[-1] != self._cells[-1]:
            return False
        if variant == "partial_top" and positions[1] != self._cells[0]:
            return False
        if variant == "partial_bottom" and positions[-2] != self._cells[-1]:
            return False
        return True

    def occurs_in(self, gc: GriddedChord) -> bool:
        max_size = len(gc)
        if max_size == 0:
            return False
        for pattern, variant in self._candidate_patterns(max_size):
            for occ in self._point_occurrences(pattern, gc):
                if self._respects_cell_constraints(gc, occ, variant):
                    return True
        return False
