from __future__ import annotations

from itertools import product
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Set, Tuple, cast

import sympy

from collections import Counter

from comb_spec_searcher.strategies.constructor.base import Constructor
from comb_spec_searcher.typing import (
    CombinatorialClassType,
    CombinatorialObjectType,
    Parameters,
    RelianceProfile,
    SubObjects,
    SubRecs,
    SubSamplers,
    SubTerms,
    Terms,
)


class ShiftedCartesianProduct(Constructor[CombinatorialClassType, CombinatorialObjectType]):
    """
    A cartesian product constructor with size shifts.

    This is used to implement Nabergall's generalized factorization (Section 2.2.4),
    where the parent class satisfies:

        T(x) = Π_i V_i(x) / x^{N_i}

    i.e. the size relation is:

        n = Σ_i (n_i - N_i).

    `shift_reliance` corresponds to the paper's set S: indices i (0-based) where the
    reliance profile for child i is truncated at n-1 (so we do not rely on size n).

    This implementation is intentionally focused on the `n` parameter; it supports
    `extra_parameters` only for the equation (variable substitution) in the same style
    as comb_spec_searcher's CartesianProduct.
    """

    def __init__(
        self,
        parent: CombinatorialClassType,
        children: Iterable[CombinatorialClassType],
        shifts: Tuple[int, ...],
        shift_reliance: Optional[Set[int]] = None,
        extra_parameters: Optional[Tuple[Dict[str, str], ...]] = None,
    ):
        self.children = tuple(children)
        if len(shifts) != len(self.children):
            raise ValueError("shifts must match number of children")
        self.shifts = tuple(int(s) for s in shifts)
        self.shift_reliance = set(shift_reliance or set())

        # Keep the same API shape as CartesianProduct for extra parameters.
        if extra_parameters is not None:
            self.extra_parameters = tuple(extra_parameters)
        else:
            self.extra_parameters = tuple(dict() for _ in self.children)

        # comb_spec_searcher expects size parameters to be nonnegative. Some classes
        # may report -1 for "no objects"; we clamp at 0 and rely on their term
        # functions returning 0 terms in that case.
        self.parent_min = max(0, cast(int, parent.minimum_size_of_object()))
        self.child_mins = tuple(
            max(0, cast(int, c.minimum_size_of_object())) for c in self.children
        )

    def can_be_equivalent(self) -> bool:  # pylint: disable=no-self-use
        """
        This constructor is not an equivalence/identity map: it changes the size
        relation by shifting, so rules using it should never be grouped as
        equivalences by the searcher.
        """
        return False

    def get_equation(
        self, lhs_func: sympy.Function, rhs_funcs: Tuple[sympy.Function, ...]
    ) -> sympy.Eq:
        res = sympy.Integer(1)
        for shift, extra_parameters, rhs_func in zip(self.shifts, self.extra_parameters, rhs_funcs):
            term = rhs_func.subs(
                {child: parent for parent, child in extra_parameters.items()},
                simultaneous=True,
            )
            if shift:
                term = term / (sympy.Symbol("x") ** shift)
            res *= term
        return sympy.Eq(lhs_func, res)

    def reliance_profile(self, n: int, **parameters: int) -> RelianceProfile:
        # We only support the size parameter `n` for now.
        profiles = []
        for idx, (child_min, shift) in enumerate(zip(self.child_mins, self.shifts)):
            hi = n + shift
            if idx in self.shift_reliance:
                hi = max(child_min, hi - 1)
            profiles.append({"n": tuple(range(child_min, hi + 1))})
        return tuple(profiles)

    def get_terms(self, parent_terms: Callable[[int], Terms], subterms: SubTerms, n: int) -> Terms:
        # `Terms` is a Counter keyed by parameter tuples (possibly empty).
        new_terms: Terms = Counter()
        ranges = []
        for idx, (child_min, shift) in enumerate(zip(self.child_mins, self.shifts)):
            hi = n + shift
            if idx in self.shift_reliance:
                hi = max(child_min, hi - 1)
            ranges.append(range(child_min, hi + 1))

        for sizes in product(*ranges):
            if sum(si - Ni for si, Ni in zip(sizes, self.shifts)) != n:
                continue

            # Multiply the children's term Counters (convolution on parameter tuples).
            prod_terms: Terms = Counter({tuple(): 1})
            for i, si in enumerate(sizes):
                child_terms = subterms[i](si)
                if not child_terms:
                    prod_terms = Counter()
                    break
                next_terms: Terms = Counter()
                for p1, v1 in prod_terms.items():
                    for p2, v2 in child_terms.items():
                        if len(p1) != len(p2):
                            raise ValueError(
                                "ShiftedCartesianProduct only supports consistent parameter dimensions."
                            )
                        next_terms[tuple(a + b for a, b in zip(p1, p2))] += v1 * v2
                prod_terms = next_terms

            new_terms += prod_terms

        return new_terms

    def get_sub_objects(
        self, subobjs: SubObjects, n: int
    ) -> Iterator[Tuple[Parameters, Tuple[List[Optional[CombinatorialObjectType]], ...]]]:
        raise NotImplementedError("Object generation is not implemented for generalized factorization.")

    def random_sample_sub_objects(
        self,
        parent_count: int,
        subsamplers: SubSamplers,
        subrecs: SubRecs,
        n: int,
        **parameters: int,
    ) -> Tuple[Optional[CombinatorialObjectType], ...]:
        raise NotImplementedError("Random sampling is not implemented for generalized factorization.")

    def equiv(
        self, other: "Constructor", data: Optional[object] = None
    ) -> Tuple[bool, Optional[object]]:
        return (
            isinstance(other, type(self))
            and self.shifts == getattr(other, "shifts", None)
            and self.shift_reliance == getattr(other, "shift_reliance", None)
            and ShiftedCartesianProduct.extra_params_equiv(
                self.extra_parameters, getattr(other, "extra_parameters", tuple())
            ),
            None,
        )

