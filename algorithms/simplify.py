"""
This module contains the SimplifyObstructionsAndRequirements class, which is used to
remove redundant obstructions and requirements and also reduce obstructions by removing factors
which must be contained.
"""

from collections import defaultdict
from itertools import product
from math import factorial
from typing import TYPE_CHECKING, Iterable

from chords import Chord, GriddedChord


def binomial(x: int, y: int) -> int:
    """Returns the binomial coefficient x choose y."""
    try:
        return factorial(x) // factorial(y) // factorial(x - y)
    except ValueError:
        return 0


class SimplifyObstructionsAndRequirements:
    """
    This class contains method for reducing and removing redundant obstructions and requirements.
    """

    def __init__(
        self,
        obstructions: tuple["GriddedChord", ...],
        requirements: tuple[tuple["GriddedChord", ...], ...],
        dimensions: tuple[int, int],
    ):
        self.obstructions = obstructions
        self.requirements = requirements
        self.dimensions = dimensions
        self.sort_obstructions()

    def remove_redundant_gridded_chords(
        self, gridded_chords: Iterable["GriddedChord"]
    ) -> tuple["GriddedChord", ...]:
        """Remove gcps that are implied by other gcps."""
        redundant_gcs = set()
        new_gridded_chords = list(gridded_chords)
        for gcp in gridded_chords:
            for gcp2 in gridded_chords:
                if gcp != gcp2 and gcp2.contains(gcp):
                    redundant_gcs.add(gcp2)
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
                tuple(req for req in req_list if req.avoids(self.obstructions))
            )
            for req_list in self.requirements
        )

    def remove_redundant_lists_requirements(self) -> None:
        """Remove requirements lists that are implied by other requirements lists."""
        indices = []
        for i, req_list_1 in enumerate(self.requirements):
            for j, req_list_2 in enumerate(self.requirements):
                if i != j and j not in indices:
                    if any(req.contains(req_list_2) for req in req_list_1):
                        indices.append(i)
        self.requirements = tuple(
            req for i, req in enumerate(self.requirements) if i not in indices
        )

    def simplify(self) -> None:
        """Simplify the obstructions and requirements using all methods until there is no change."""
        curr_obs = None
        curr_reqs = None
        while curr_obs != self.obstructions or curr_reqs != self.requirements:
            curr_obs = self.obstructions
            curr_reqs = self.requirements
            self.simplify_once()
            self.sort_requirements()
            self.sort_obstructions()

    def simplify_once(self) -> None:
        """Do one pass of all the different simplify methods."""
        self.remove_redundant_obstructions()
        self.remove_redundant_requirements()
        self.remove_redundant_lists_requirements()
        self.remove_factors_from_obstructions()

    def sort_requirements(self) -> None:
        """Orders the requirements and removes duplicates."""
        self.requirements = tuple(
            sorted(set(tuple(sorted(set(req_list))) for req_list in self.requirements))
        )

    def sort_obstructions(self) -> None:
        """Orders the obstructions and removes duplicates."""
        self.obstructions = tuple(sorted(set(self.obstructions)))

    def remove_factors_from_obstructions(self) -> None:
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
        cells = ob.find_active_cells()
        for factor in ob.find_factors(self.point_rows()):
            if self.implied_by_requirements(factor):
                cells.difference_update(factor.find_active_cells())
        return ob.sub_gridded_cayley_perm(cells)

    def point_rows(self) -> set[int]:
        """
        Returns the point rows of the tiling.

        # TODO: be passed from the tiling in the init to avoid duplicated code?
        """
        point_rows = set()
        counter_dict: dict[int, int] = defaultdict(int)
        for ob in self.obstructions:
            if ob.pattern in (CayleyPermutation([0, 1]), CayleyPermutation([1, 0])):
                if ob.positions[0][1] == ob.positions[1][1]:
                    counter_dict[ob.positions[0][1]] += 1
        for row, count in counter_dict.items():
            n = len(self.cells_in_row(row))
            if 2 * binomial(n, 2) + 2 * n == count:
                point_rows.add(row)
        return point_rows

    def cells_in_row(self, row: int) -> set[tuple[int, int]]:
        """Returns the set of active cells in the given row."""
        cells = set()
        for cell in self.active_cells():
            if cell[1] == row:
                cells.add(cell)
        return cells

    def active_cells(self) -> set[tuple[int, int]]:
        """Returns the set of active cells in the tiling.
        (Cells are active if they do not contain a point obstruction.)"""
        active_cells = set(
            product(range(self.dimensions[0]), range(self.dimensions[1]))
        )
        for ob in self.obstructions:
            if len(ob) == 1:
                active_cells.discard(ob.positions[0])
        return active_cells

    def implied_by_requirement(
        self, gcp: "GriddedCayleyPerm", req_list: Iterable["GriddedCayleyPerm"]
    ) -> bool:
        """Check whether a gridded Cayley permutation is implied by a requirement."""
        return all(req.contains_gridded_cperm(gcp) for req in req_list)

    def implied_by_requirements(self, gcp: "GriddedCayleyPerm") -> bool:
        """Check whether a gridded Cayley permutation is implied by the requirements."""
        return any(
            self.implied_by_requirement(gcp, req_list) for req_list in self.requirements
        )