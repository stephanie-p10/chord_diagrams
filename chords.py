from comb_spec_searcher import CombinatorialClass, CombinatorialObject
from tilings import Tiling
from tilings import GriddedPerm
from typing import Callable, Dict, FrozenSet, Iterable, Iterator, List, Optional, Tuple, Union
from itertools import combinations, islice, tee, product, chain
from permuta.misc import UnionFind
import json

DIR_EAST = 0
DIR_NORTH = 1
DIR_WEST = 2
DIR_SOUTH = 3
DIR_NONE = -1
DIRS = [DIR_EAST, DIR_NORTH, DIR_WEST, DIR_SOUTH]

Cell = Tuple[int, int]
Position = Tuple[Cell, ...]

class ChordException(Exception):
    def __init__(self, msg=None):
        if msg == None:
            msg == "Chord error"
        else:
            msg == "Chord error: " + msg
        super().__init__(msg)

class Chord(Tuple):
    def __new__(cls, iterable: Iterable[int] = ()) -> "Chord":
        """Creates and validates a new Chord instance. 
        Examples:
        >>> Chord((0, 0, 1, 1))
        (0, 0, 1, 1)

        >>> Chord.from_iterable_validated((0,0,0))
        Chord has uneven length

        >>> Chord.from_iterable_validated((0,0,0,0))
        Duplicate chord 0

        >>> Chord.from_iterable_validated((1,1,0,0))
        Incorrect indexing, chord 1 came before chord 0

        >>> Chord.from_iterable_validated((0, None))
        Traceback (most recent call last):
        ...
        TypeError: 'None' object is not an integer
        """

        # Ensure length is even (uses bitwise operator)
        if len(iterable) % 2 ==1:
            raise ChordException(f"Chord has uneven length: {len(iterable)}")
        
        # Loops through and checks conditions for a valid chord
        max_chord_seen = -1
        unpaired_chords = set()
        for num in iterable:
            if not isinstance(num, int):
                raise TypeError(f"'{num}' object is not an integer")
            # The chord num skiped a number from the biggest chord seen so far.
            if num > max_chord_seen + 1:
                raise ChordException(f"Incorrect indexing, chord {num} came before chord {num-1}")
            # The chord num is the next chord after the biggest chord seen so far.
            elif num == max_chord_seen + 1:
                max_chord_seen = num
                # First time seeing num, so it is currently unseen
                unpaired_chords.add(max_chord_seen) 
            elif num <= max_chord_seen:
                if num in unpaired_chords:
                    # num was seen once, now we are seeing it a second time, so the chord has become paired.
                    unpaired_chords.remove(num) 
                else:
                    raise ChordException(f"Duplicate chord {num}")
        if len(unpaired_chords) != 0:
            raise ChordException(f"The chord {unpaired_chords.pop()} is unpaired.")

        return tuple.__new__(cls, iterable)

    def __init__(self, iterable: Iterable[int] = ()) -> None:
        # Cache for data used when finding occurrences of self in a perm
        self._cached_pattern_details: Optional[List[Tuple[int, int, int, int]]] = None
        self._length = len(iterable)//2
        self._pattern = iterable

        # Creates dictionary of chord number and vertices it connects
        # (0, 1, 2, 1, 3, 3, 2, 0) becomes {0: (0, 7), 1: (1, 3), 2: (2, 6), (3: (4, 5)}
        self._chord = {}
        for i,c in enumerate(iterable):
            if c not in self._chord:
                self._chord[c] = (i, -1)
            else:
                self._chord[c] = (self._chord[c][0],i)

    def __len__(self) -> int:
        return self._length
    
    # s added
    def get_pattern(self) -> list[int]:
        """Returns the pattern of the chord as a list"""
        return self._pattern
    
    # sCP
    def concat(self, *others: "Chord") -> "Chord":
        """Return the concatenation of two or more Chords in the given order.

        Examples:
            >>> Chord((0, 0, 1, 1)).concat(Chord((0, 1, 1, 0)))
            Chord((0, 0, 1, 1, 2, 3, 3, 2))
        """
        result = list(self)
        for other in others:
            result.extend(val + self._length for val in other)
        return Chord(result)

    # sCP
    def remove_chord(self, selected: Optional[int] = None) -> "Chord":
        """Return the Chord acquired by removing i-th chord (0-indexed). It
        defaults to the last chord.

        Examples:
            >>> Chord((0, 1, 0, 1, 2, 2)).remove_chord()
            Chord((0, 1, 0, 1)))
            >>> Chord((0, 1, 1, 0, 2, 2)).remove_chord(2)
            Chord((0, 1, 1, 0)))
            >>> Chord((0, 1, 1, 0, 2, 3, 3, 2)).remove_chord(3)
            Chord((0, 1, 1, 0, 2, 2)))
        """
        if selected is None:
            if self._length == 0:
                return self
            selected = self._length - 1
        assert 0 <= selected and selected < self._length
        return Chord(
            list(val if val < selected else val - 1 for val in self if val != selected)
        )

    def remove_index(self, index: Optional[int] = None) -> "Chord":
        """Return the Chord acquired by removing a chord at a specified index. It
        defaults to the index of the diagram

        Examples:
            >>> Chord((0, 1, 0, 1, 2, 2)).remove_index()
            Chord((0, 1, 0, 1)))
            >>> Chord((0, 1, 1, 0, 2, 2)).remove_index(2)
            Chord((0, 0, 1, 1)))
            >>> Chord((0, 1, 1, 0, 2, 3, 3, 2)).remove_index(3)
            Chord((0, 0, 1, 2, 2, 1)))
        """
        selected = -1
        if index is None:
            if self._length == 0:
                return self
            selected = self[-1] # defaults to the last index
        else:
            selected = self[index]
        return self.remove_chord(selected)
    
    # sCP, sAsk: it feels like there is a better way to do this
    # sToDo: test with this function commented out, then get rid
    def standardize(self, original_list: list) -> "Chord":
        """Return a Chord based on the original list ordering.
        Assumes valid chord is inputted

        Examples:
            >>> Chord(()).standardize([-1, 1, -1, 2, 2, 1])
            Chord((0, 1, 0, 2, 2, 1))
            >>> Chord(()).standardize([0, 1, 0, 2, 2, 1, 5, 5])
            Chord((0, 1, 0, 2, 2, 1, 3, 3))
        """
        unique_elements = []
        for num in original_list:
            if num not in unique_elements:
                unique_elements.append(num)
        
        unique_elements.sort()
    
        result_list = []
        for num in original_list:
            result_list.append(unique_elements.index(num))
        return Chord(result_list)
    
    @classmethod
    def to_standard(cls, iterable):
        # Keeps track of standard order chords items in iterable will be mapped to.
        chord_map = {}
        # Keeps track of what chord the next new item in iterable will be mapped to.
        next_index = 0
        
        # Maps each item in iterable to its appropriate chord.
        for item in iterable:
            if item not in chord_map:
                chord_map[item] = next_index
                next_index += 1
        
        # Puts iterable in standard order using chord_map to map items to their chord values.
        return cls(tuple(chord_map[item] for item in iterable))

    # sCP: similar method for this was implemented in occurances_in
    def _contains(self, patt: "Chord") -> bool:
        """Determines if a chord contains an instance of patt
        Examples:
            >>> Chord((0, 1, 1, 0)).contains(Chord((0, 0)))
            True
            >>> Chord((0, 1, 1, 0)).contains(Chord((0, 1, 0, 1)))
            False"""
        if isinstance(patt, Chord):
            return len(patt.occurrences_in(self)) > 0
        raise TypeError("patt must be a Chord")
    
    # sCN: changed function in structure (just calls occurrences_in)
    def contains(self, *chords: "Chord") -> bool:
        """
        Check if self contains all patts.
        Examples:
            >>> Chord((0,0,1,1)).contains(Chord((0,0)))
            True
            >>> Chord((0,0,1,1)).contains(Chord((0,0,1,1)), Chord((0,1,0,1)))
            False
        """
        return all(self._contains(patt) for patt in chords)

    # sCN: changed function majorly in structure
    def occurrences_in(self, patt: "Chord", colour_self: Iterator = None, 
                       colour_patt: Iterator[int] = None) -> Iterator[Tuple[int, ...]]:
        """Find all indices of occurrences of self in patt. 
        Examples:
            >>> list(Chord((0, 1, 0, 1)).occurrences_in(Chord((0,1,2,0,1,2))))
            [(0, 1), (1, 2), (0, 2)]
            >>> list(Chord((0,0)).occurrences_in(Chord((0, 1, 1, 0))))
            [(0,), (1,)]
            >>> list(Chord().occurrences_in(Chord((0, 1, 1, 0))))
            [()]
        """
        self_has_colour = (colour_self != None)
        patt_has_colour = (colour_patt != None)

        if (not self_has_colour and patt_has_colour) or (not patt_has_colour and self_has_colour):
            return [()]
        
        if self_has_colour and len(colour_self) != len(self) * 2:
            raise Exception("Incorrect colour length for self")
        
        if patt_has_colour and len(colour_patt) != len(patt) * 2:
            raise Exception("Incorrect colour length for pattern given")

        if len(self) == 0 or len(self) > len(patt):
            return [()]
        
        instances = []

        for subchord in combinations(patt._chord, len(self)):
            subchord_indices_in_patt = []
            for i in subchord:
                subchord_indices_in_patt.extend(patt._chord[i])
            subchord_indices_in_patt.sort()

            subchord_patt = []
            for i in subchord_indices_in_patt:
                subchord_patt.append(patt._pattern[i])

            if Chord.to_standard(subchord_patt) == self:
                append = True
                if patt_has_colour and self_has_colour:
                    for i in range(len(subchord_indices_in_patt)):
                        if colour_patt[subchord_indices_in_patt[i]] != colour_self[i]:
                            append = False
                if append:
                    instances.append(subchord)

        return instances

    # sCP, added examples
    def avoids(self, *chords: "Chord") -> bool:
        """
        Check if self avoids all chords.
        Examples:
        >>> Chord((0, 0, 1, 1)).avoids(Chord((0, 0, 1, 1)), Chord((0, 1, 0, 1))))
        False
        >>> Chord((0, 1, 1, 0, 2, 3, 3, 2)).avoids(Chord((0, 1, 0, 1)), Chord((0, 1, 2, 2, 1, 0))))
        True
        >>> Chord((0, 1, 1, 2, 0, 3, 3, 2)).avoids(Chord((0, 1, 1, 0))))
        False
        """
        return all(not self._contains(patt) for patt in chords)

    # sCP
    def insert(self, m: Optional[int] = None, n: Optional[int] = None) -> "Chord":
        """Return the Chord acquired by inserting a new chord at a specified index. The 
        new chord is inserted at m-points from the root and n-points from the root. 
        If m=0, then the new chord's left side is the root. The default value is another chord inserted
        at the end.

        Examples:
            >>> Chord((0, 0, 1, 1)).insert()
            Chord((0, 0, 1, 1, 2, 2))
            >>> Chord((0, 1, 0, 1)).insert(0, 0)
            Chord((0, 0, 1, 2, 1, 2))
            >>> Chord((0, 1, 2, 2, 0, 1)).insert(0, 4)
            Chord((0, 1, 2, 3, 3, 0, 1, 2))
        """
        listlen = 2*self._length
        if m is None:
            return self.concat(Chord((0, 0)))
        if n is None:
            n = listlen
        
        assert (0 <= m) and (m <= n) and (n <= listlen)
        result = list(self)
        result.insert(m, -1)
        result.insert(n+1, -1)
        return Chord.to_standard(result)

    # sCN: added default, which has the functionality of get_chord, returning whole chord
    def get_chords(self, selected: list = None) -> "Chord":
        """Return the Chord based on the selected list of Chords. If no argument 
        is given, it will retrun the whole chord.
        Examples:
            >>> Chord((0, 0, 1, 1)).get_chords([0])
            Chord((0, 0))
            >>> Chord((0, 1, 0, 1)).get_chords([0, 1])
            Chord((0, 1, 0, 1))
            >>> Chord((0, 1, 1, 0, 2, 2)).get_chords()
            Chord((0, 1, 1, 0, 2, 2)))
            >>> Chord((0, 1, 1, 0, 2, 2)).get_chords([1, 2])
            Chord((0, 0, 1, 1)))
        """
        if selected == None:
            return self
        
        return Chord.to_standard([val for val in self if val in selected])

    # sCP, sCN: changed name from chord -> chord_dict, added example
    @property
    def chord_dict(self) -> dict[(int)]:
        """Returns the chord dictionary
        Example:
            >>> Chord((0, 1, 1, 0, 2, 2)).chord_dict
            {0: (0, 3), 1: (1, 2), 2: (4, 5)}"""
        return self._chord
    
    # sToDo: all of the following up to connected (these methods will probably be fun)
    def crossing(self):
        pass
    def nesting(self):
        pass
    def k_crossing(self):
        pass
    def connected(self):
        pass



# sToDo: GriddedChord should initialize with a chord instead of a list, and this should be fixed in methods as well
class GriddedChord(CombinatorialObject):
    def __init__(
        self, chord: Chord = Chord(()), positions: Iterable[Cell] = ()
    ) -> None:
        self._chord = chord
        self._patt = chord.get_pattern()
        self._pos = tuple(positions)
        if len(self._patt) != len(self._pos):
            raise ValueError("Pattern and position list have unequal lengths.")
        self._cells: FrozenSet[Cell] = frozenset(self._pos)

    @classmethod
    def single_cell(cls, chord: Chord, cell: Cell) -> "GriddedChord":
        """Construct a gridded chord where the chords are all located in a
        single cell.

        Examples:
        >>> GriddedChord.single_cell(Chord((0, 0, 1, 1, 2, 2)), (0,0)) 
        GriddedChord(Chord((0, 0, 1, 1, 2, 2)), ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)))
        >>> GriddedChord.single_cell(Chord((0, 1, 0, 2, 2, 1)), (1, 1))
        GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1) ))
        """
        return cls(chord, (cell for _ in range(int(len(chord.get_pattern())))))

    @classmethod
    def empty_chord(cls) -> "GriddedChord":
        """Construct the empty gridded chord.

        Examples:
        GriddedChord.empty_chord() == GriddedChord()
        GriddedChord.empty_chord() == GriddedChord(Chord(), ())"""
        return cls(Chord(), ())

    @classmethod 
    # sCN: point_chord -> single_chord
    def single_chord(cls, cells: Iterable[Cell]) -> "GriddedChord":
        """Construct the single gridded chord using the cells given. If only one cell 
        is given, will construct the single gridded chord with both ends in the same cell
        
        Examples:
        >>> GriddedChord.single_chord((0,0))
        GriddedChord((0, 0), ((0, 0), (0, 0)))
        >>> GriddedChord.single_chord((1, 0), (1, 1))
        GriddedChord((0, 0), ((1, 0), (1, 1)))"""
        if len(cells) == 1:
            return cls(Chord((0,0)), (cells[0],cells[0]))
        elif len(cells) == 2:
            if cells[0][0] != cells[1][0]:
                raise Exception("inconsistant columns given")
            return cls(Chord((0,0)), (cells[0],cells[1]))
        raise ValueError("incorrect number of cells given")

    def occupies(self, cell: Cell) -> bool:
        """Checks if the gridded chord has a point in the given cell.
        
        Examples:
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))).occupies((1, 1))
        True
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))).occupies((0, 1))
        False
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (0,1), (1,2))).occupies((0, 0))
        True
        """
        return cell in self._cells

    # sToDo: could be generalized to recognize patterns that are not indexed the same
    def occurrences_in(self, other: "GriddedChord") -> Iterator[Tuple[int, ...]]:
        """Returns all occurrences of self in other.

        Examples:
        gc1 = GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (0,1), (1,2)))
        gc2 = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), ((0, 0), (1, 0), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4)))
        gc3 = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), ((0, 0), (1, 1), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 2), (4, 4)))
        >>> gc1.occurrences_in(gc2)
        (0, 2)
        >>> gc1.occurrences_in(gc3)
        (0, 1), (0, 2) 
        """
        # uses method from Chord class, with patterns given (so the occurrences have to match in positions)
        yield from self._chord.occurrences_in(other._chord, self._pos, other._pos)

    def occurs_in(self, other: "GriddedChord") -> bool:
        """Checks if self occurs in other.
        Examples:
        gc1 = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), 
                                ((0, 0), (1, 0), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4)))
        gc2 = GriddedChord(Chord((0, 1, 0, 2, 2, 1)), 
                                ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1)))
        gc3 = GriddedChord(Chord((0, 1, 0, 1)),
                                ((0,0), (1,1), (0,1), (1,2)))
        gc4 = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), 
                                ((0, 0), (1, 1), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 2), (4, 4)))

        >>> gc3.occurs_in(gc1)
        True
        >>> gc3.occurs_in(gc2)
        False
        >>> gc3.occurs_in(gc4)
        True
        >>> gc1.occurs_in(gc3)
        False
        """
        return any(self.occurrences_in(other))

    def contains(self, *patts: "GriddedChord") -> bool:
        """Return true if self contains an occurrence of any of patts.
        Examples:
        gc = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), 
                               ((0, 0), (1, 0), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4)))
        >>> gc.contaions(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3))))
        True
        >>> gc.contains(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 3), (0, 1)))
        False
        >>> gc.contains(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3))), 
                        GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 3), (0, 1)))
        True
        """
        return any(any(True for _ in patt.occurrences_in(self)) for patt in patts)
    
    def avoids(self, *patts: "GriddedChord") -> bool:
        """Return true if self avoids all of the patts.
        
        Examples:
        gc = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), 
                               ((0, 0), (1, 0), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4)))
        >>> gc.contaions(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3))))
        True
        >>> gc.contains(GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 3), (0, 1)))
        False
        >>> gc.contains(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3))), 
                        GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 3), (0, 1)))
        False
        """
        return not self.contains(*patts)

    #sCN: the old method was perm specific
    def contradictory(self) -> bool:
        """Checks if the points of the griddedchord contradict the chord.

        Checks if for every 0 <= i < j < n, 
            if i and j are on the same chord, their x positions are the same
            if one chord comes before the other, the x position is less or equal
            the y positions are in increasing order 

        Examples:
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (0,1), (1,2))).contradictory()
        False
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (1,0), (1,2))).contradictory()
        True
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((1,0), (0,1), (1,1), (0,2))).contradictory()
        True
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((1, 0), (0, 1), (1, 1), (1, 2))).contradictory()
        True
        """
        """Old method (for perms, does not work for chords):
        return any(
            (l_x > r_x or (l_v < r_v and l_y > r_y) or (l_v > r_v and l_y < r_y)) or (l_v == r_v and l_x != r_x)
            for idx, (j_chord, (j_x, j_y)) in enumerate(self)
            for i_chord, (i_x, i_y) in islice(self, idx)
        )"""
        # loops over combinations of chords with i < j in self
        for idx, (chordj, (xj, yj)) in enumerate(self):
            for chordi, (xi, yi) in islice(self, idx):
                # positions for points on the same chord must be in the same column
                if (chordj == chordi) and (xj != xi):
                    return True
                
                # The x positions of increasing chords must increase (or stay the same)
                if ((chordj < chordi) and (xj > xi)) or ((chordj > chordi) and (xj < xi)):
                    return True
                    
                # y positions must be in increasing order of the chord
                if (yj < yi):
                    return True
                    
        return False
        
    # sCN: old method was perm specific
    def remove_cells(self, cells: Iterable[Cell]) -> "GriddedChord":
        """Remove any chords in the cell given and return a new gridded
        chord.
        
        Examples:
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))).remove_cells((1, 1))
        GriddedChord()
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 1), (0, 1), (1, 2)).remove_cells((0, 0))
        GriddedChord(Chord((0, 0)), ((1, 1), (1, 2)))
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 1), (0, 1), (1, 2)).remove_cells((0, 0), (1, 1))
        GriddedChord()"""
        cells = set(cells)
        chords_to_remove = [chord for chord, pos in self if pos in cells]
        remaining_chords = []
        remaining_positions = []

        for chord, pos in self: 
            # add the chords and positions that are not being removed to lists
            if chord not in chords_to_remove:
                remaining_chords.append(chord)
                remaining_positions.append(pos)
        
        #print(Chord.to_standard(remaining_chords))

        return type(self)(
            Chord.to_standard(remaining_chords),
            remaining_positions
        )

    def points_in_cell(self, cell: Cell) -> Iterator[int]:
        """Yields the indices of the points in the cell given.
        
        Examples: 
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))).points_in_cell()
        (0, 1, 2, 3, 4, 5)
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 1), (1, 1), (0, 1), (1, 1), (1, 1), (1, 1))).points_in_cell()
        (2, 3, 4, 5)"""
        return (i for i, (_, pos) in enumerate(self) if pos == cell)

    def isolated_cells(self) -> Iterator[Cell]:
        """Yields the cells that contain only one chord of the gridded
        chord and are in their own row and column (besides the other end of the chord).

        Examples:
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).isolated_cells()
        (0, 0)
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (2, 2), (2, 2), (1, 3))).isolated_cells()
        (0, 0), (1, 1), (2, 2)
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 2), (2, 2), (1, 2))).isolated_cells()
        (0, 0)
        """
        return set(
            (x1, y1)
            for i, (chord1, (x1, y1)) in enumerate(self)
            if not any(
                x1 == x2 or y1 == y2 for j, (chord2, (x2, y2)) in enumerate(self) 
                if (i != j and chord1 != chord2))
        )

    # sCN: added check for chords being in same column
    def is_isolated(self, indices: Iterable[int]) -> bool:
        """Checks if the cells at the indices do not share a row or column with
        any other cell in the gridded chord.
        
        Examples:
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([0])
        True
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([0, 1])
        False
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (1, 1), (1, 1), (1, 1))).is_isolated([1])
        False
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 0), (1, 1), (0, 0), (2, 2), (2, 2), (1, 2))).is_isolated([3])
        False"""
        indices = set(indices)
        return not any(
            x_j == x_i or y_j == y_i
            for j, (chord_j, (x_j, y_j)) in enumerate(self)
            for i, (chord_i, (x_i, y_i)) in enumerate(self)
            if i not in indices and j in indices and chord_i != chord_j
        )

    # sCN: returns most extreme chord in each case, not point
    def forced_point_index(self, cell: Cell, direction: int) -> int:
        """Search in the cell given for the chord with the strongest force with
        respect to the given force.
        
        Examples:
        """
        if self.occupies(cell):
            indices = self.points_in_cell(cell)
            if direction == DIR_EAST:
                return max(chord for chord, _ in self)
            if direction == DIR_NORTH:
                return self._patt[max(indices)]
            if direction == DIR_WEST:
                return min(chord for chord, _ in self)
            if direction == DIR_SOUTH:
                return self._patt[max(indices)]
            raise ValueError("You're lost, no valid direction")
        raise ValueError("The gridded chord does not occupy the cell")

    # sAsk: not sure what this is doing, why it is important, and why I should care. Moving on... 
    def forced_point_of_requirement(
        self, gps: Tuple["GriddedChord", ...], indices: Tuple[int, ...], direction: int
    ) -> Optional[Tuple[int, int]]:
        """
        Return the pair (x, y) where x is the gridded chord in gps that is
        farthest in the given direction in self, and y is index of the forced point
        with respect to the gps and indices. If gps is avoided, then
        return None.
        """

        def directionmost(i1: int, i2: int) -> int:
            """return the directionmost between i1 and i2."""
            if direction == DIR_EAST:
                return i1 if self._patt[i1] > self._patt[i2] else i2
            if direction == DIR_NORTH:
                return max(i1, i2)
            if direction == DIR_WEST:
                return i1 if self._patt[i1] < self._patt[i2] else i2
            if direction == DIR_SOUTH:
                return min(i1, i2)
            raise ValueError("You're lost, no valid direction")

        res: Optional[Tuple[int, int]] = None
        for idx, gp in enumerate(gps):
            forced_index_in_patt = indices[idx]
            for occurrence in gp.occurrences_in(self):
                if res is None:
                    res = idx, occurrence[forced_index_in_patt]
                else:
                    new_res = directionmost(res[1], occurrence[forced_index_in_patt])
                    if res[1] != new_res:
                        res = idx, new_res
        return res

    def get_points_col(self, col: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord in the column col."""
        return ((val, i) for i, (val, (x, _)) in enumerate(self) if x == col)

    def get_points_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord in the row."""
        return ((val, i) for i, (val, (_, y)) in enumerate(self) if y == row)

    def get_points_below_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord below the row."""
        return ((val, i) for i, (val, (_, y)) in enumerate(self) if y < row)

    def get_points_above_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord above the row."""
        return ((val, i) for i, (val, (_, y)) in enumerate(self) if y > row)

    def get_points_left_col(self, col) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord left of column col."""
        return ((val, i) for i, (val, (x, _)) in enumerate(self) if x < col)
    
    def get_points_right_col(self, col: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord right of column col."""
        return ((val, i) for i, (val, (x, _)) in enumerate(self) if x > col)

    def get_subchord_left_col(self, col: int) -> "GriddedChord":
        """Returns the gridded subchord of points left of column col."""
        gen1, gen2 = tee((i for chord, i in self.get_points_left_col(col)), 2)
        return type(self)(
            Chord.to_standard([self._patt[i] for i in gen1])._pattern, 
            (self._pos[i] for i in gen2)
        )

    # sCN: used to be from indices, now from chords (previously called get_gridded_chord_at_indices(self, indices: Iterable[int]))
    def get_subgrid_at_chords(self, chords: Iterable[int]) -> "GriddedChord":
        """
        Returns the subgridded chord that contains only the point at the given
        indices. Indices must be sorted.
        """
        indices = [idx for idx, (chord, pos) in enumerate(self) if chord in chords]
        gen1, gen2 = tee(indices, 2)
        return type(self)(
            Chord.to_standard([self._patt[i] for i in gen1])._pattern,
            (self._pos[i] for i in gen2),
        )

    # sCN: used cleaner list notation
    def remove_chord(self, index) -> "GriddedChord":
        avoid_chord = self._patt[index]
        new_chords = [chord for chord, _ in self if chord != avoid_chord]
        new_positions = [pos for chord, pos in self if chord != avoid_chord]

        new_patt = Chord.to_standard(new_chords)._pattern

        return GriddedChord(new_patt, new_positions)
            
    # sCN: gets all chords with at least one endpoint in cells - should this be non inclusive?
    def get_gridded_chord_in_cells(self, cells: Iterable[Cell]) -> "GriddedChord":
        """Returns the subgridded chord with chords with endpoints in cells."""
        chords_to_keep = [chord for chord, pos in self if pos in cells]
        patt = []
        positions = []

        for chord, pos in self:
            if chord in chords_to_keep:
                patt.append(chord)
                positions.append(pos)

        new_grid = GriddedChord(Chord.to_standard(patt)._pattern, positions)


        return new_grid

    # sCN: previous method returned the restrictions on the values for a point inserted 
    # into a given cell for permutations. Changed to accomadate chords.
    def get_bounding_box(self, row1: int, row2: int) -> Tuple[int, int, int, int]:
        """Determines the range possible indices for a chord inserted with sink in row1 and source in row2
        Returns in order: (min chord, max chord, min index source, max index source, min index sink, max index sink)"""
        assert(row1 <= row2)

        rows_below_row1 = list(self.get_points_below_row(row1))
        rows_above_row1 = list(self.get_points_above_row(row1))
        # finds minimun index a source of a new chord in row1 can have, by checking what the largest chord below row1 is
        min_index_row1 = max([idx for chord, idx in rows_below_row1] + [-1]) + 1
        # finds maximum index a source of a new chord in row1, by checking what the smallest chord above row1 is.
        max_index_row1 = min([idx for chord, idx in rows_above_row1] + [len(self) * 2])

        rows_below_row2 = list(self.get_points_below_row(row2))
        rows_above_row2 = list(self.get_points_above_row(row2))
        # finds minimum index a sink of a new chord in row2 can have, checking what the largest chord below row2 is, 
        # adding 1 since we assume the source of the chord will be established below the sink.
        min_index_row2 = max([idx + 1 for chord, idx in rows_below_row2] + [min_index_row1]) + 1
        # finds maximum index a new sink chord can have in row2, checking what the smallest chord above row2 is,
        # adding 1 since we assume the source of the chord will be established below the sink.
        max_index_row2 = min([idx for chord, idx in rows_above_row2] + [max_index_row1]) + 1

        return (min_index_row1, max_index_row1, min_index_row2, max_index_row2)

    # sCN: insert_point changed to insert_chord, functionality changed appropriately
    def insert_chord(self, col: int, row1: int, row2: int) -> Iterator["GriddedChord"]:
        """Insert a new point into cell of the gridded chord, such that the
        point is added to the underlying pattern with the position at the cell.
        Yields all gridded chords where the point has been mixed into the points
        in the cell."""
        min_index_source, max_index_source, min_index_sink, max_index_sink = self.get_bounding_box(col, row1, row2)
        # sAsk: There is probably a better way to do this
        for source in range(min_index_source, max_index_source + 1):
            for sink in range(min(source, min_index_sink), max_index_sink):
                yield self.insert_specific_chord(col, row1, row2, source, sink)

    # sCN: insert_specific_point -> insert_specific_chord, changed numebr of parameters, adapted to chords appropriately
    #sToDo: code check for vaild gridded chord required
    def insert_specific_chord(self, col: int, row1: int, row2: int, source: int, sink: int) -> "GriddedChord":
        """Insert the chord chord in column col, with sink in row1 at index sink, and source in row2 at index source"""
        patt = self._patt
        patt = patt.insert(source, sink)
        positions = list(self._pos)
        positions.insert(source, (col, row1))
        positions.insert(sink, (col, row2))
        return GriddedChord(patt._pattern, positions)
    
    #sCN: changed to work with chords
    def remove_chord(self, chord: int) -> "GriddedChord":
        """Remove the point at index from the gridded chord."""
        patt = Chord.to_standard([item for item, pos in self if item != chord])._pattern
        pos = [pos for item, pos in self if item != chord]
        return type(self)(patt, pos)

    def all_subchords(self, proper: bool = True) -> Iterator["GriddedChord"]:
        """Yields all gridded subchords."""
        # loops through all sizes of subchords, taking into account if we are looking for proper subchords or not
        for subchord_length in range(len(self) if proper else len(self) + 1):
            for subchords in combinations(range(len(self)), subchord_length):
                yield type(self)(
                    Chord.to_standard([chord for chord, pos in self if chord in subchords])._pattern,
                    (pos for chord, pos in self if chord in subchords)
                )

    # this method does not seem relavent, also don't know what it is doing
    '''
    def extend(self, c: int, r: int) -> Iterator["GriddedChord"]:
        """Add n+1 to all possible positions in chord and all allowed positions given
        that placement."""
        n = len(self)
        if n == 0:
            yield from (
                GriddedChord((0,0), ((x, y),)) for x in range(c) for y in range(r)
            )
        else:
            min_y = max(y for _, (_, y) in self)
            yield from (
                self.insert_specific_chord((x, y), index, n)
                for index in range(n + 1)
                for x in range(
                    self._pos[index - 1][0] if index > 0 else 0,
                    self._pos[index][0] + 1 if index < n else c,
                )
                for y in range(min_y, r)
            )'''

    # this looks like it works? not sure what to test on it, if it is necessary, so it is being left as is
    def apply_map(self, cell_mapping: Callable[[Cell], Cell]) -> "GriddedChord":
        """Map the coordinates to a new list of coordinates according to the
        cell_mapping given."""
        return type(self)(self._patt._pattern, [cell_mapping(cell) for cell in self._pos])
    
    # sCN: is_point_chord -> is_single_chord
    def is_single_chord(self) -> bool:
        """Checks if the gridded chord is of length 1."""
        return len(self) == 1

    # sToDo: why is this function just calling another? redundant?
    def is_localized(self) -> bool:
        """Check if the gridded chord occupies only a single cell."""
        return self.is_single_cell()

    def is_single_cell(self) -> bool:
        """Check if the gridded chord occupies only a single cell."""
        return len(self._cells) == 1

    def is_single_row(self) -> bool:
        """Check if the gridded chord occupies only a single row."""
        return len(set(y for (_, y) in self._cells)) == 1
    
    def is_single_col(self) -> bool:
        """Check if the gridded chord occupies only a single column"""
        return len(set(x for (x, _) in self._cells)) == 1

    def is_empty(self) -> bool:
        """Check if the gridded chord is the gridded chord."""
        return len(self) == 0

    def is_interleaving(self) -> bool:
        """Check if the gridded chord occupies two cells that are in the
        same row or column."""
        xs: Dict[int, int] = {}
        ys: Dict[int, int] = {}
        for x, y in self._cells:
            # get keys x from xs, y from ys, if not found return y, x resp.. Is return value not y or x?
            if xs.get(x, y) != y or ys.get(y, x) != x:
                return True
            # add to keys
            xs[x] = y
            ys[y] = x
        return False
    
    # sToDo: not sure what this is doing, need to look more at union find algorithm
    # sAsk: what is unionfind
    def factors(self) -> List["GriddedChord"]:
        """Return a list containing the factors of a gridded chord.
        A factor is a sub gridded chord that is isolated on its own rows
        and columns."""
        uf = UnionFind(len(self._pos))
        for j, (_, (x_r, y_r)) in enumerate(self):
            for i, (_, (x_l, y_l)) in enumerate(islice(self, j)):
                if x_l == x_r or y_l == y_r:
                    uf.unite(i, j)
        # Collect the connected factors of the cells
        all_factors: Dict[int, List[Cell]] = {}
        for i, cell in enumerate(self.pos):
            x = uf.find(i)
            if x in all_factors:
                all_factors[x].append(cell)
            else:
                all_factors[x] = [cell]
        factor_cells = list(set(cells) for cells in all_factors.values())
        return [self.get_gridded_chord_in_cells(comp) for comp in factor_cells]

    def compress(self) -> List[int]:
        """Compresses the gridded chord into a list of integers.
        It starts with a list of the values in the chord. The rest is
        the list of positions flattened."""
        return list(chain(self._patt, chain.from_iterable(self._pos)))

    @classmethod
    def decompress(cls, arr: Iterable) -> "GriddedChord":
        """Decompresses a list of integers in the form outputted by the
        compress method and constructs an Obstruction."""
        length = len(arr) // 3
        iterator = iter(arr)
        return cls(
            tuple(next(iterator) for _ in range(length)), ((next(iterator), next(iterator)) for _ in range(length))
        )

    # looks about right, have not tested.
    def apply_chord_map_to_cell(
        self, chord_mapping: Callable[[Chord], Chord], cell: Cell
    ) -> "GriddedChord":
        """Apply a chord map to the subchord within a cell."""
        subchord = [val for val, pos in self if pos == cell]
        st = Chord.to_standard(subchord)
        back_map = dict(zip(st, subchord))
        new_subchord = (back_map[val] for val in chord_mapping(st))
        return type(self)(
            (next(new_subchord) if pos == cell else val for val, pos in self), self.pos
        )

    def to_jsonable(self) -> dict:
        """Returns a dictionary object which is JSON serializable representing
        a GriddedChord."""
        return {"patt": self._patt, "pos": self._pos}

    @classmethod
    def from_json(cls, jsonstr: str) -> "GriddedChord":
        """Returns a GriddedChord object from JSON string."""
        jsondict = json.loads(jsonstr)
        return cls.from_dict(jsondict)

    @classmethod
    def from_dict(cls, jsondict: dict) -> "GriddedChord":
        """Returns a GriddedChord object from a dictionary loaded from a JSON
        serialized GriddedChord object."""
        return cls(jsondict["patt"], map(tuple, jsondict["pos"]))  # type: ignore

    @property
    def patt(self) -> Chord:
        return self._patt

    @property
    def pos(self) -> Tuple[Cell, ...]:
        return self._pos
    
    def ascii_plot(self) -> str:
        max_x = max(cell[0] for cell in self._pos)
        max_y = max(cell[1] for cell in self._pos)
        res = ""

        def points_in_col(i):
            return sum(1 for cell in self._pos if cell[0] == i)//2

        def points_in_row(j):
            return sum(1 for cell in self._pos if cell[1] == j)

        row_boundary = (
            "+" + "+".join("-" * points_in_col(i) for i in range(max_x + 1)) + "+"
        )
        col_boundary = (
            "|" + "|".join(" " * points_in_col(i) for i in range(max_x + 1)) + "|"
        )

        # Starting at the highest y position, until -1 is reached, step by -1
        # creates empty grid
        for j in range(max_y, -1, -1):
            res += (
                "\n".join(
                    [row_boundary] + [col_boundary for i in range(points_in_row(j))]
                )
                + "\n"
            )
        res += row_boundary
        
        for idx, (val, (x, y)) in enumerate(self):
            # insert into this spot:
            # (val + x + 1) is the horizontal index. val is chords to left, and
            #               x + 1 counts number of col boundaries to left
            # (len(self) + max_y) - (idx + y) is vertical index.
            #               idx is points below, and y + 1 counts number of - below
            insert = (val + x + 1) + ((len(self)*2 + max_y) - (idx + y)) * (
                len(row_boundary) + 1)
            res = res[:insert] + "●" + res[insert + 1 :]

        #for chord, (x, y) in self:
            

        return res
    
    def __len__(self) -> int:
        return len(self._patt)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({tuple(self._patt)!r}, {self._pos})"

    def __str__(self) -> str:
        return f"{self._patt}: {', '.join(str(c) for c in self._pos)}"

    def __hash__(self) -> int:
        return hash(self._patt) ^ hash(self._pos)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return False
        return self._patt == other._patt and self._pos == other._pos

    def __lt__(self, other: "GriddedChord") -> bool:
        return (self._patt, self._pos) < (other._patt, other._pos)

    def __contains__(self, other: "GriddedChord") -> bool:
        return next((True for _ in other.occurrences_in(self)), False)

    def __iter__(self) -> Iterator[Tuple[int, Cell]]:
        return zip(self._patt, self._pos)

