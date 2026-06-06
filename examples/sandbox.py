import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from src.common.chords import Chord, GriddedChord
from src.common.tiling import Tiling
from itertools import product

listA = [1, 2, 3]
listB = [4, 5]
listC = [6]
req_list = [listA, listB, listC]

print(3 in listA)
