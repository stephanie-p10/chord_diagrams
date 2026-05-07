from itertools import combinations
from math import factorial
from typing import Iterable, Iterator, List, Sequence, Set, Tuple, TypeVar

T = TypeVar("T")


def intersection_reduce(sets: Iterable[Iterable[T]]) -> Set[T]:
    """Return the intersection of an iterable of iterables.

    Returns a list (stable order is not guaranteed).
    """
    it = iter(sets)
    try:
        first = set(next(it))
    except StopIteration:
        return set()
    for s in it:
        first.intersection_update(s)
    return first


def union_reduce(sets: Iterable[Iterable[T]]) -> Set[T]:
    """Return the union of an iterable of iterables.

    Returns a list (stable order is not guaranteed).
    """
    res = set()
    for s in sets:
        res.update(s)
    return res


def multinomial(*parts: int) -> int:
    """Compute multinomial coefficient for given part sizes."""
    n = sum(parts)
    denom = 1
    for p in parts:
        denom *= factorial(p)
    return factorial(n) // denom


def partitions_iterator(items: Sequence[T]) -> Iterator[List[List[T]]]:
    """Yield all set partitions of a sequence as lists of blocks.

    This is a simple recursive implementation intended for small inputs.
    """
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for partition in partitions_iterator(rest):
        # put first in its own block
        yield [[first]] + partition
        # insert first into each existing block
        for i in range(len(partition)):
            new_part = [block[:] for block in partition]
            new_part[i].append(first)
            yield new_part

