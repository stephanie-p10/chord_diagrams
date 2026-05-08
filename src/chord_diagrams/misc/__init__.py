
"""Small shared constants and helpers.

Currently this module primarily defines direction constants used by placement
and navigation code:

- `DIR_EAST`, `DIR_NORTH`, `DIR_WEST`, `DIR_SOUTH`: the four cardinal directions
- `DIR_NONE`: sentinel for “no direction”
- `DIRS`: the ordered list of the four cardinal directions
"""

DIR_EAST = 0
DIR_NORTH = 1
DIR_WEST = 2
DIR_SOUTH = 3
DIR_NONE = -1
DIRS = [DIR_EAST, DIR_NORTH, DIR_WEST, DIR_SOUTH]

__all__ = [
    "DIRS",
    "DIR_EAST",
    "DIR_NORTH",
    "DIR_WEST",
    "DIR_SOUTH",
    "DIR_NONE",
]