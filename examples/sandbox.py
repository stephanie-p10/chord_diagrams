'''import sys
from pathlib import Path
 
_src_root = Path(__file__).resolve().parents[2]  # .../src
if str(_src_root) not in sys.path:
    sys.path.insert(0, str(_src_root))

from src.chord_diagrams.misc import DIR_EAST, DIR_SOUTH

from src.chord_diagrams.chords import Chord, GriddedChord
from src.chord_diagrams.tiling import Tiling'''
from itertools import product

listA = [1, 2, 3]
listB = [4, 5]
listC = [6]
req_list = [listA, listB, listC]

print(3 in listA)
