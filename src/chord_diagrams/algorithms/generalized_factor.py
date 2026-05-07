import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from collections import defaultdict
from itertools import chain, combinations
from typing import TYPE_CHECKING, Dict, Iterable, Iterator, List, Optional, Set, Tuple

from permuta.misc import UnionFind
from chords import GriddedChord, Chord
from tilings.assumptions import ComponentAssumption, TrackingAssumption
from tilings.misc import partitions_iterator

if TYPE_CHECKING:
    from tilings import Tiling


Cell = Tuple[int, int]
ReqList = Tuple[GriddedChord, ...]

__all__ = ("GeneralizedFactor",)


class GeneralizedFactor:
    """
    Implements a restricted version of Nabergall's generalized factorization (Section 2.2.4).

    The full paper definition works with *factorizable covers* Q of the nonempty cells,
    where the induced subtilings may overlap on intersection subtilings U^(ij) of size 0 or 1,
    and the resulting size relation is shifted by the overlap sizes.

    In this codebase we currently only implement the most basic special case:
    - all overlaps are on a single common, atom subtiling U (|U| = 1),
    - every factor is of the form A_i = C_i ∪ U where C_i are the strong components
      (row/col + obstruction/requirement/linkage) of the complement.

    """
    def __init__(self, tiling: "Tiling") -> None:
        self._tiling = tiling
        self._active_cells = tiling.active_cells
        nrow = tiling.dimensions[1]
        ncol = tiling.dimensions[0]
        self._cell_unionfind = UnionFind(nrow * ncol)
        self._components: Optional[Tuple[Set[Cell], ...]] = None
        self._factors_obs_and_reqs: Optional[
            List[
                Tuple[
                    Tuple[GriddedChord, ...],
                    Tuple[ReqList, ...],
                    Tuple[TrackingAssumption, ...],
                ]
            ]
        ] = None

    def _cell_to_int(self, cell: Cell) -> int:
        nrow = self._tiling.dimensions[1]
        return cell[0] * nrow + cell[1]

    def _int_to_cell(self, i: int) -> Cell:
        nrow = self._tiling.dimensions[1]
        return (i // nrow, i % nrow)

    def _get_cell_representative(self, cell: Cell) -> Cell:
        """
        Return the representative of a cell in the union find.
        """
        i = self._cell_to_int(cell)
        return self._cell_unionfind.find(i)  # type: ignore

    def _unite_cells(self, cells: Iterable[Cell]) -> None:
        """
        Put all the cells of `cells` in the same component of the UnionFind.
        """
        cell_iterator = iter(cells)
        try:
            c1 = next(cell_iterator)
        except StopIteration:
            return
        c1_int = self._cell_to_int(c1)
        for c2 in cell_iterator:
            c2_int = self._cell_to_int(c2)
            self._cell_unionfind.unite(c1_int, c2_int)

    def _unite_assumptions(self) -> None:
        """
        For each TrackingAssumption unite all the positions of the gridded perms.
        """
        for assumption in self._tiling.assumptions:
            if isinstance(assumption, ComponentAssumption):
                for cells in assumption.cell_decomposition(self._tiling):
                    self._unite_cells(cells)
            else:
                for gp in assumption.gps:
                    self._unite_cells(gp.pos)

    def _unite_obstructions(self) -> None:
        """
        For each obstruction unite all the position of the obstruction.
        """
        for ob in self._tiling.obstructions:
            self._unite_cells(ob.pos)

    def _unite_requirements(self) -> None:
        """
        For each requirement unite all the cell in all the requirements of the
        list.
        """
        for req_list in self._tiling.requirements:
            req_cells = chain.from_iterable(req.pos for req in req_list)
            self._unite_cells(req_cells)

    @staticmethod
    def _same_row_or_col(cell1: Cell, cell2: Cell) -> bool:
        """
        Test if `cell1` and `cell2` or in the same row or columns
        """
        return cell1[0] == cell2[0] or cell1[1] == cell2[1]

    def _unite_rows_and_cols(self) -> None:
        """
        Unite all the active cell that are on the same row or column.
        """
        cell_pair_to_unite = (
            c
            for c in combinations(self._active_cells, r=2)
            if self._same_row_or_col(c[0], c[1])
        )
        for c1, c2 in cell_pair_to_unite:
            self._unite_cells((c1, c2))

    def _unite_all(self) -> None:
        """
        Unite all the cells that share an obstruction, a requirement list,
        a row or a column.
        """
        self._unite_obstructions()
        self._unite_requirements()
        self._unite_assumptions()
        self._unite_rows_and_cols()

    def get_components(self) -> Tuple[Set[Cell], ...]:
        """
        Returns the tuple of all the components. Each component is set of
        cells.
        """
        if self._components is not None:
            return self._components
        self._unite_all()
        all_components: Dict[Cell, Set[Cell]] = defaultdict(set)
        for cell in self._active_cells:
            rep = self._get_cell_representative(cell)
            all_components[rep].add(cell)
        self._components = tuple(all_components.values())
        return self._components

    def _get_factors_obs_and_reqs(
        self,
    ) -> List[
        Tuple[
            Tuple[GriddedChord, ...], Tuple[ReqList, ...], Tuple[TrackingAssumption, ...]
        ],
    ]:
        """
        Returns a list of all the irreducible factors of the tiling.
        Each factor is a tuple (obstructions, requirements)
        """
        if self._factors_obs_and_reqs is not None:
            return self._factors_obs_and_reqs
        if self._tiling.is_empty():
            return [((GriddedChord(Chord(), []),), tuple(), tuple())]
        
        factors = []
        for component in self.get_components():
            #print(component)
            obstructions = tuple(
                ob for ob in self._tiling.obstructions if ob.pos[0] in component
            )
            requirements = tuple(
                req for req in self._tiling.requirements if req[0].pos[0] in component
            )
            # TODO: consider skew/sum assumptions
            assumptions = tuple(
                ass.__class__(gp for gp in ass.gps if gp.pos[0] in component)
                for ass in self._tiling.assumptions
            )
            #print(obstructions, requirements)
            factors.append(
                (
                    obstructions,
                    requirements,
                    tuple(set(ass for ass in assumptions if ass.gps)),
                )
            )
        self._factors_obs_and_reqs = factors
        return self._factors_obs_and_reqs

    def factorable(self) -> bool:
        """
        Returns `True` if the tiling has more than one factor.
        """
        return len(self.get_components()) > 1

    # ------------------------
    # Generalized factorization
    # ------------------------

    def _cells_in_common_patterns(self) -> Set[Cell]:
        """
        Return cells participating in any obstruction/requirement/linkage.

        Used to identify candidate overlap subtile cells U that are weakly isolated
        from the rest of the tiling (i.e. no pattern spans U and its complement).
        """
        cells: Set[Cell] = set()
        for ob in self._tiling.obstructions:
            cells.update(ob.pos)
        for req_list in self._tiling.requirements:
            for req in req_list:
                cells.update(req.pos)
        for linkage in self._tiling.linkages:
            cells.update(linkage)
        for ass in self._tiling.assumptions:
            for gp in getattr(ass, "gps", ()):
                cells.update(gp.pos)
        return cells

    def _is_weakly_non_interacting(self, s1: Set[Cell], s2: Set[Cell]) -> bool:
        """
        Check the paper's *weakly non-interacting* condition for s1 and s2:
        no obstruction, requirement list, linkage, or assumption involves cells from both.
        """
        if not s1 or not s2:
            return True

        # obstructions
        for ob in self._tiling.obstructions:
            pos = set(ob.pos)
            if pos & s1 and pos & s2:
                return False

        # requirement lists
        for req_list in self._tiling.requirements:
            cells = set(chain.from_iterable(req.pos for req in req_list))
            if cells & s1 and cells & s2:
                return False

        # linkages
        for linkage in self._tiling.linkages:
            cells = set(linkage)
            if cells & s1 and cells & s2:
                return False

        # assumptions (treated like patterns on their gp positions)
        for assumption in self._tiling.assumptions:
            if isinstance(assumption, ComponentAssumption):
                for cells in assumption.cell_decomposition(self._tiling):
                    cells = set(cells)
                    if cells & s1 and cells & s2:
                        return False
            else:
                for gp in getattr(assumption, "gps", ()):
                    cells = set(gp.pos)
                    if cells & s1 and cells & s2:
                        return False
        return True

    @staticmethod
    def _is_strongly_non_interacting(s1: Set[Cell], s2: Set[Cell]) -> bool:
        if not s1 or not s2:
            return True
        for c1 in s1:
            for c2 in s2:
                if c1[0] == c2[0] or c1[1] == c2[1]:
                    return False
        return True

    def generalized_factorizations(
        self,
    ) -> Iterator[Tuple[Tuple[Tuple[Cell, ...], ...], Tuple[int, ...], Set[int]]]:
        """
        Yield generalized factorizations as (parts, shifts, shift_reliance).

        - parts is a cover Q of active cells.
        - shifts is the tuple N from the paper, one per part.
        - shift_reliance is the paper's set S, represented as 0-based child indices.
        """
        active = set(self._active_cells)
        if len(active) <= 1:
            return

        # Candidate U: a subset of active cells that is weakly non-interacting with its complement
        # and whose induced subtiling is an atom (|U| = 1).
        # We only try single-cell U first, then (rarely) pairs.
        candidates: List[Set[Cell]] = [{c} for c in active]
        if len(active) <= 12:
            candidates += [set(p) for p in combinations(active, 2)]

        for U in candidates:
            rest = active - U
            if not rest:
                continue
            if not self._is_weakly_non_interacting(U, rest):
                continue

            U_tiling = self._tiling.sub_tiling(tuple(sorted(U)))
            if not U_tiling.is_atom():
                continue

            # Split the remainder by strong interaction (row/col) *plus* weak interactions.
            # We approximate this using the existing union-find components on the full tiling,
            # then refine by removing U.
            comps = [set(comp) - U for comp in self.get_components()]
            comps = [c for c in comps if c]
            if len(comps) <= 1:
                continue

            # Build cover parts A_i = C_i ∪ U
            parts = [tuple(sorted(c | U)) for c in comps]

            # Ensure each part has something outside U.
            if any(len(set(p) - U) == 0 for p in parts):
                continue

            # Compute shifts N as in the paper for the case where every intersection is U.
            nU = U_tiling.minimum_size_of_object()
            shifts = [0] + [nU for _ in range(1, len(parts))]

            # Compute shift reliance S using minimum sizes (sufficient condition for |V_{N_i}|=0).
            children = [self._tiling.sub_tiling(p) for p in parts]
            shift_reliance: Set[int] = set()
            for i in range(len(children)):
                for j in range(len(children)):
                    if i == j:
                        continue
                    if children[j].minimum_size_of_object() > shifts[j]:
                        shift_reliance.add(i)
                        break

            yield tuple(parts), tuple(shifts), shift_reliance
            return

    def factors(self) -> Tuple["Tiling", ...]:
        """
        Returns all the irreducible factors of the tiling.
        """
        #print(self._get_factors_obs_and_reqs())
        return tuple(
            self._tiling.__class__(
                obstructions=f[0],
                requirements=f[1],
                assumptions=tuple(sorted(f[2])),
                simplify=False,
            )
            for f in self._get_factors_obs_and_reqs()
        )

    def reducible_factorizations(self) -> Iterator[Tuple["Tiling", ...]]:
        """
        Iterator over all reducible factorization that can be obtained by
        grouping of irreducible factors.

        Each factorization is a list of tiling.

        For example if T = T1 x T2 x T3 then (T1 x T3) x T2 is a possible
        reducible factorization.
        """
        min_comp = self._get_factors_obs_and_reqs()
        for partition in partitions_iterator(min_comp):
            factors = []
            for part in partition:
                obstructions, requirements, assumptions = zip(*part)
                factors.append(
                    self._tiling.__class__(
                        obstructions=chain(*obstructions),
                        requirements=chain(*requirements),
                        assumptions=tuple(sorted(chain(*assumptions))),
                        simplify=False,
                    )
                )
            yield tuple(factors)

