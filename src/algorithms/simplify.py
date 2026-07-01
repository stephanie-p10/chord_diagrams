"""
This module contains the SimplifyObstructionsAndRequirements class, which is used to
remove redundant obstructions and requirements and also reduce obstructions by removing factors
which must be contained.
"""

from collections import defaultdict
from itertools import combinations, product
from math import factorial
from typing import TYPE_CHECKING, Dict, FrozenSet, Iterable, List, Optional, Set, Tuple

import sys
from pathlib import Path

_src_root = Path(__file__).resolve().parents[2]  # .../src
if str(_src_root) not in sys.path:
    sys.path.insert(0, str(_src_root))

from src.common.chords import Chord, GriddedChord

Cell = Tuple[int, int]
Linkage = Tuple[Cell, ...]


def binomial(x: int, y: int) -> int:
    """Returns the binomial coefficient x choose y."""
    try:
        return factorial(x) // factorial(y) // factorial(x - y)
    except ValueError:
        return 0


def _linkage_as_frozenset(link: Iterable[Cell]) -> FrozenSet[Cell]:
    return frozenset(link)


def _normalize_linkages(
    linkages: Iterable[Iterable[Cell]],
) -> Tuple[FrozenSet[Cell], ...]:
    """Return linkages as deduplicated frozensets in stable order."""
    seen: Set[FrozenSet[Cell]] = set()
    result: List[FrozenSet[Cell]] = []
    for link in linkages:
        fs = _linkage_as_frozenset(link)
        if fs and fs not in seen:
            seen.add(fs)
            result.append(fs)
    return tuple(result)


def _linkage_hypergraph_connected(linkages: Tuple[FrozenSet[Cell], ...]) -> bool:
    """True when linkages share cells so their intersection graph is connected."""
    if len(linkages) <= 1:
        return True
    visited = {0}
    queue = [0]
    while queue:
        idx = queue.pop()
        for j, other in enumerate(linkages):
            if j in visited:
                continue
            if linkages[idx] & other:
                visited.add(j)
                queue.append(j)
    return len(visited) == len(linkages)


def _linkage_rows_of(link: FrozenSet[Cell]) -> Set[int]:
    return {cell[1] for cell in link}


def _linkage_row_partition_connected(linkages: Tuple[FrozenSet[Cell], ...]) -> bool:
    """Condition 2: each row meets at most one linkage; row groups are connected."""
    if len(linkages) <= 1:
        return True

    row_owner: dict[int, int] = {}
    for idx, link in enumerate(linkages):
        for cell in link:
            row = cell[1]
            if row in row_owner and row_owner[row] != idx:
                return False
            row_owner[row] = idx

    rows_by_link = [_linkage_rows_of(link) for link in linkages]
    visited = {0}
    queue = [0]
    while queue:
        idx = queue.pop()
        rows_i = rows_by_link[idx]
        if not rows_i:
            continue
        min_i, max_i = min(rows_i), max(rows_i)
        for j, rows_j in enumerate(rows_by_link):
            if j in visited or not rows_j:
                continue
            min_j, max_j = min(rows_j), max(rows_j)
            if max_i + 1 >= min_j and max_j + 1 >= min_i:
                visited.add(j)
                queue.append(j)
    return len(visited) == len(linkages)


def _linkage_union_covers(
    target: FrozenSet[Cell],
    candidates: Tuple[FrozenSet[Cell], ...],
) -> Optional[Tuple[FrozenSet[Cell], ...]]:
    """Return a subset of candidates whose union equals target, if one exists."""
    contained = tuple(link for link in candidates if link and link <= target)
    if not contained:
        return None
    for size in range(1, len(contained) + 1):
        for subset in combinations(contained, size):
            if set().union(*subset) == target:
                return subset
    return None


def _linkage_apriori_finite_connected(
    obstructions: Tuple[GriddedChord, ...],
    requirements: Tuple[Tuple[GriddedChord, ...], ...],
    linkages: Tuple[FrozenSet[Cell], ...],
    linkage: FrozenSet[Cell],
    active_cells: Iterable[Cell],
    dimensions: Tuple[int, int],
) -> bool:
    """Condition 3: on T \\ {L}, only finitely many L-restrictions and all connected."""
    # Singleton linkages are kept as seeds for strategy-induced expansion
    # (e.g. requirement placement stretching linkages to new cells).
    if len(linkage) <= 1:
        return False

    remaining = tuple(
        tuple(sorted(link)) for link in linkages if link != linkage
    )
    from src.common.tiling import Tiling

    tiling_without = Tiling(
        obstructions,
        requirements,
        remaining,
        simplify=False,
        expand=False,
        remove_empty_rows_and_cols=False,
        derive_empty=False,
    )

    cells = tuple(sorted(linkage))
    max_size = min(2 * len(cells) + 2, 12)
    subdiagrams: Set[Tuple[Tuple[int, ...], Tuple[Cell, ...]]] = set()

    for size in range(0, max_size + 1):
        for chord in Chord.of_length(size):
            for gc in GriddedChord.all_grids(chord, active_cells):
                if not tiling_without.contains(gc):
                    continue
                if not any(cell in linkage for cell in gc.pos):
                    continue
                sub = gc.get_subchord_in_cells(list(cells))
                subdiagrams.add((sub.patt, sub.pos))
                if not sub.is_connected(list(cells)):
                    return False

    return len(subdiagrams) > 0


def _linkage_deletable_as_union(
    linkage: FrozenSet[Cell],
    linkages: Tuple[FrozenSet[Cell], ...],
) -> bool:
    """Conditions 1 and 2: linkage is the union of other connected linkages."""
    others = tuple(link for link in linkages if link != linkage)
    cover = _linkage_union_covers(linkage, others)
    if cover is None or len(cover) < 2:
        return False
    return _linkage_hypergraph_connected(cover) or _linkage_row_partition_connected(cover)


def _linkage_can_delete(
    linkage: FrozenSet[Cell],
    linkages: Tuple[FrozenSet[Cell], ...],
    obstructions: Tuple[GriddedChord, ...],
    requirements: Tuple[Tuple[GriddedChord, ...], ...],
    active_cells: Iterable[Cell],
    dimensions: Tuple[int, int],
) -> bool:
    """Return True if linkage satisfies a sufficient deletion condition."""
    if not linkage:
        return True

    if _linkage_deletable_as_union(linkage, linkages):
        return True

    others = tuple(link for link in linkages if link != linkage)

    if any(linkage < other for other in others):
        return True

    # Keep singleton linkages: they carry connectivity intent for later expansion.
    if len(linkage) == 1:
        return False

    return _linkage_apriori_finite_connected(
        obstructions,
        requirements,
        linkages,
        linkage,
        active_cells,
        dimensions,
    )


def _remove_redundant_linkages(
    linkages: Iterable[Iterable[Cell]],
    obstructions: Tuple[GriddedChord, ...],
    requirements: Tuple[Tuple[GriddedChord, ...], ...],
    active_cells: Iterable[Cell],
    dimensions: Tuple[int, int],
) -> Tuple[Linkage, ...]:
    """Remove linkages that are redundant under Nabergall's sufficient conditions."""
    normalized = _normalize_linkages(linkages)
    changed = True
    while changed:
        changed = False
        for link in normalized:
            if _linkage_deletable_as_union(link, normalized):
                normalized = tuple(l for l in normalized if l != link)
                changed = True
                break
        if changed:
            continue
        for link in normalized:
            if _linkage_can_delete(
                link,
                normalized,
                obstructions,
                requirements,
                active_cells,
                dimensions,
            ):
                normalized = tuple(l for l in normalized if l != link)
                changed = True
                break
    return tuple(tuple(sorted(link)) for link in normalized)


# s TODO: a lot of this class has redundant code compared to tilings
class SimplifyObstructionsAndRequirements:
    """
    This class contains method for reducing and removing redundant obstructions and requirements.
    """

    def __init__(
        self,
        obstructions: Tuple["GriddedChord", ...],
        requirements: Tuple[Tuple["GriddedChord", ...], ...],
        dimensions: Tuple[int, int],
        active_cells: Tuple[Tuple[int, int]],
        empty_cells: Tuple[Tuple[int, int]],
        linkages: Tuple[Tuple[Tuple[int, int], ...], ...] = tuple(),
    ):
        self.obstructions = obstructions
        self.requirements = requirements
        self.dimensions = dimensions
        self.active_cells = active_cells
        self.empty_cells = empty_cells
        self.linkages = linkages
        self.sort_obstructions()

    def remove_redundant_gridded_chords(
        self, gridded_chords: Iterable["GriddedChord"]
    ) -> Tuple["GriddedChord", ...]:
        """Remove gcs that are implied by other gcs."""
        redundant_gcs = set()
        new_gridded_chords = list(gridded_chords)
        for gc in gridded_chords:
            for gc2 in gridded_chords:
                if gc != gc2 and gc2.contains(gc):
                    redundant_gcs.add(gc2)
        for gcs in redundant_gcs:
            new_gridded_chords.remove(gcs)
        return tuple(new_gridded_chords)

    def remove_redundant_obstructions(self) -> None:
        """Remove obstructions that are implied by other obstructions."""
        self.obstructions = self.remove_redundant_gridded_chords(self.obstructions)

    def remove_redundant_requirements(self) -> None:
        """Remove requirements that are implied by other requirements in the same list."""
        
        self.requirements = tuple(
            self.remove_redundant_gridded_chords(
                tuple(req for req in req_list if req.avoids(*self.obstructions))
            )
            for req_list in self.requirements
        )

    def remove_redundant_lists_requirements(self) -> None:
        """Remove requirements lists that are implied by other requirements lists."""
        indices = []
        for i, req_list_1 in enumerate(self.requirements):
            if self.requirement_implied_by_some_requirement(
                req_list_1,
                [
                    req_list_2
                    for j, req_list_2 in enumerate(self.requirements)
                    if i != j and j not in indices
                ],
            ):
                indices.append(i)
        self.requirements = tuple(
            req for i, req in enumerate(self.requirements) if i not in indices
        )
    
    def remove_redundant_cells(self) -> None:
        """
        Here we simplify a cell to have a point-obstruction if it 
        originally has a chord obstruction in it and all other active cells 
        in the row have a chord obstruction between them and the cell.
        """

        active_cells = set(self.active_cells)
        if not active_cells:
            return

        obs_set = set(self.obstructions)

        def has_single_chord_obstruction_in_cells(c1: Tuple[int, int], c2: Tuple[int, int]) -> bool:
            # A single-chord obstruction is represented by pattern (0,0) placed in two cells.
            # The cell could have the source or sink of the chord so check both.
            return (
                GriddedChord(Chord((0, 0)), (c1, c2)) in obs_set
                or GriddedChord(Chord((0, 0)), (c2, c1)) in obs_set
            )

        cells_to_remove: Set[Tuple[int, int]] = set()

        # removes cells that are empty due to chord obstructions in the row
        for cell in active_cells:
            # Must obstruct a chord entirely within the cell.
            if not has_single_chord_obstruction_in_cells(cell, cell):
                continue
            row_cells = {c for c in active_cells if c[1] == cell[1] and c != cell}
            if not row_cells:
                continue
            # Must obstruct a chord between this cell and every other active cell in the row.
            if all(has_single_chord_obstruction_in_cells(cell, other) for other in row_cells):
                cells_to_remove.add(cell)

        # removes cells that are empty due to RGF properties of chord diagrams.
        required_cells = [cell for cell in active_cells if self.implied_by_requirements(GriddedChord(Chord((0,)), (cell,)))]
        for cell in active_cells:
            southwest_active_cells = [c2 for c2 in active_cells if (c2[0] <= cell[0] and c2[1] <= cell[1]) and not 
                                                                   (c2[0] == cell[0] and c2[1] == cell[1])]
            southeast_required_cells = [c2 for c2 in required_cells if (c2[0] > cell[0] and c2[1] < cell[1])]
            if len(southwest_active_cells) == 0 and len(southeast_required_cells) > 0:
                cells_to_remove.add(cell)

        if not cells_to_remove:
            return

        # Add a point obstruction for each removed cell
        point_obs = tuple(
            GriddedChord.single_cell(Chord((0,)), cell) for cell in sorted(cells_to_remove)
        )
        self.obstructions = tuple(sorted(set(self.obstructions + point_obs)))

        self.update_cells_status()
        #print("updated obstructions", self.obstructions[0])

    def requirement_implied_by_some_requirement(
        self,
        requirement: Tuple["GriddedChord", ...],
        requirements: Iterable[Tuple["GriddedChord", ...]],
    ) -> bool:
        """Check if one of the requirements implies the containment of requirement."""
        return any(
            self.requirement_implied_by_requirement(requirement, req)
            for req in requirements
        )
    
    @staticmethod
    def requirement_implied_by_requirement(
        requirement: Tuple["GriddedChord", ...],
        other_requirement: Tuple["GriddedChord", ...],
    ) -> bool:
        """Check if the containment of other implies containment of requirement."""
        return all(
            any(other_gc.contains(gc) for gc in requirement)
            for other_gc in other_requirement
        )

    def remove_redundant_linkages(self) -> None:
        """Remove linkages implied by other linkages or a priori connectivity."""
        if not self.linkages:
            return
        self.linkages = _remove_redundant_linkages(
            self.linkages,
            self.obstructions,
            self.requirements,
            self.active_cells,
            self.dimensions,
        )

    def simplify(self) -> None:
        """Simplify the obstructions and requirements using all methods until there is no change."""
        curr_obs = None
        curr_reqs = None
        curr_links = None
        while (
            curr_obs != self.obstructions
            or curr_reqs != self.requirements
            or curr_links != self.linkages
        ):
            
            #print(self.obstructions[0])
            curr_obs = self.obstructions
            curr_reqs = self.requirements
            curr_links = self.linkages
            self.simplify_once()
            #print("after simplify once", self.obstructions[0])
            self.sort_requirements()
            #print("after sort requirements", self.obstructions[0])
            self.sort_obstructions()
            #print("after sort obs", self.obstructions[0])
            #print()

    def simplify_once(self) -> None:
        """Do one pass of all the different simplify methods."""
        self.remove_redundant_linkages()
        self.remove_redundant_cells()
        #print("before removing redundant obs", self.obstructions[0])
        self.remove_redundant_obstructions()
        #print("after removing redundant obs", self.obstructions[0])
        #print("after remove redundant", self.obstructions)
        self.remove_redundant_requirements()
        #print("after removing requirements", self.obstructions[0])
        self.remove_redundant_lists_requirements()
        #print("after removing requirements lists", self.obstructions[0])
        # maybe add this back in, not sure what it is doing, but it skrews with obstructions in tilings tests
        #self.remove_factors_from_obstructions() 
        #print("after remove factors", self.obstructions)

    def sort_requirements(self) -> None:
        """Orders the requirements and removes duplicates."""
        self.requirements = tuple(
            sorted(set(tuple(sorted(set(req_list))) for req_list in self.requirements))
        )

    def sort_obstructions(self) -> None:
        """Orders the obstructions and removes duplicates."""
        self.obstructions = tuple(sorted(set(self.obstructions)))

    '''def remove_factors_from_obstructions(self) -> None:
        """Removes factors from all of the obstructions."""
        self.obstructions = tuple(
            self.remove_factors_from_obstruction(ob) for ob in self.obstructions
        )

    def remove_factors_from_obstruction(
        self, ob: "GriddedChord"
    ) -> "GriddedChord":
        """
        Removes factors from a single obstruction:

        Splits an obstruction into its factors and removes the factors that are
        implied by the requirements.
        """
        cells = set(ob.get_active_cells())
        for factor in ob.find_factors(self.point_rows()):
            if self.implied_by_requirements(factor):
                cells.difference_update(factor.get_active_cells())
        return ob.get_subchord_in_cells(cells)

    # sToDo: fix where this comes from: it is duplicated from tilings?
    def point_rows(self) -> Set[int]:
        """
        Returns the point rows of the tiling. -> when only one value is in the row

        #TODO: be passed from the tiling in the init to avoid duplicated code?
        # s TODO this is wrong, point cells need to be calculated differently with expand
        """
        point_rows = set()
        counter_dict: Dict[int, int] = defaultdict(int)
        for ob in self.obstructions:
            if ob.patt in (Chord([0, 1]), Chord([1, 0])):
                if ob.pos[0][1] == ob.pos[1][1]:
                    counter_dict[ob.pos[0][1]] += 1
        for row, count in counter_dict.items():
            n = len(self.cells_in_row(row))
            if 2 * binomial(n, 2) + 2 * n == count:
                point_rows.add(row)
        return point_rows

    def cells_in_row(self, row: int) -> Set[Tuple[int, int]]:  # this looks like it should be in the tilings repo
        """Returns the set of active cells in the given row."""
        cells = set()
        for cell in self.active_cells():
            if cell[1] == row:
                cells.add(cell)
        return cells'''

    def update_cells_status(self):  # should be in tilings repo
        """Returns the set of active cells in the tiling.
        (Cells are active if they do not contain a point obstruction.)"""
        active_cells = set(self.active_cells)
        for ob in self.obstructions:
            # Empty cells are marked by a point obstruction `Chord((0,))`.
            if ob.is_point() and ob.pos:
                active_cells.discard(ob.pos[0])
        self.active_cells = list(active_cells)
        self.empty_cells = tuple(
                cell
                for cell in product(range(self.dimensions[0]), range(self.dimensions[1]))
                if cell not in active_cells
            ) 
        

    def implied_by_requirement(
        self, gc: "GriddedChord", req_list: Iterable["GriddedChord"]
    ) -> bool:
        """Check whether a gridded chord diagram is implied by a requirement."""
        return all(req.contains(gc) for req in req_list)

    def implied_by_requirements(self, gc: "GriddedChord") -> bool:
        """Check whether a gridded chord diagram is implied by the requirements."""
        return any(
            self.implied_by_requirement(gc, req_list) for req_list in self.requirements
        )
    

    