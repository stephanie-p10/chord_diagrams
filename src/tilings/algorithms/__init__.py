"""Compatibility re-exports for code that expects `tilings.algorithms`."""

from chord_diagrams.algorithms.factor import Factor
from chord_diagrams.algorithms.obstruction_inferral import SubobstructionInferral


# These variants are used only by some strategies; keep as aliases for now.
FactorWithInterleaving = Factor
FactorWithMonotoneInterleaving = Factor

__all__ = [
    "Factor",
    "FactorWithInterleaving",
    "FactorWithMonotoneInterleaving",
    "SubobstructionInferral",
]

