from typing import Iterable, Iterator, Tuple


class TrackingAssumption:
    """Minimal stand-in for `tilings.assumptions.TrackingAssumption`."""

    def __init__(self, gps: Iterable) -> None:
        self.gps = tuple(gps)


class ComponentAssumption(TrackingAssumption):
    """A tracking assumption that can decompose into cell components."""

    def cell_decomposition(self, tiling) -> Iterator[Tuple[Tuple[int, int], ...]]:
        # Fallback: treat each gp position-set as its own component.
        for gp in self.gps:
            pos = getattr(gp, "pos", ())
            if pos:
                yield tuple(pos)


class SumComponentAssumption(ComponentAssumption):
    pass


class SkewComponentAssumption(ComponentAssumption):
    pass


__all__ = [
    "TrackingAssumption",
    "ComponentAssumption",
    "SumComponentAssumption",
    "SkewComponentAssumption",
]

