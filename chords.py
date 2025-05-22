from comb_spec_searcher import CombinatorialClass, CombinatorialObject
from tilings import Tiling
from tilings import GriddedPerm
from typing import Callable, Dict, FrozenSet, Iterable, Iterator, List, Optional, Tuple, Union
from itertools import combinations, islice, tee

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
        if len(iterable) & 1:
            raise ChordException("Chord has uneven length")
        
        # Loops through and checks conditions for a valid chord
        max_chord_seen = -1
        unpaired_chords = set()
        for num in iterable:
            if not isinstance(num, int):
                raise TypeError(f"'{num}' object is not an integer")
            if num > max_chord_seen + 1:
                raise ChordException(f"Incorrect indexing, chord {num} came before chord {num-1}")
            elif num == max_chord_seen + 1:
                max_chord_seen = num
                # since this is the first time seeing num, it is currently an unpaird chord.
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
        # (1,2,1,3,3,2) -> {(1, (1, 3)), (2, (2, 6)), (3, (4, 5))}
        self._chord = {}
        for i,c in enumerate(iterable):
            if c not in self._chord:
                self._chord[c] = (i, -1)
            else:
                self._chord[c] = (self._chord[c][0],i)

    def __len__(self) -> int:
        return self._length
    
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

    # sToDo: fix examples
    def remove_index(self, index: Optional[int] = None) -> "Chord":
        """Return the Chord acquired by removing a chord at a specified index. It
        defaults to the index of the diagram

        Examples:
            >>> Chord((0, 0, 1, 1)).remove_index()
            Chord((0, 0))
            >>> Chord((0, 0, 1, 1, 2, 2)).remove_index(0)
            Chord((0, 0, 1, 1))
            >>> Chord((0, 1, 1, 0, 2, 3, 3, 2)).remove_index(1)
            Chord((0, 0, 1, 2, 2, 1))
        """
        selected = -1
        if index is None:
            if self._length == 0:
                return self
            selected = self[-1] # defaults to the last index
        else:
            selected = self[index]
        return self.remove_chord(selected)
    
    # sToDo: fix examples; sCP
    def remove_chord(self, selected: Optional[int] = None) -> "Chord":
        """Return the Chord acquired by removing i-th chord (0-indexed). It
        defaults to the last chord.

        Examples:
            >>> Chord((0, 0, 1, 1)).remove_chord()
            Chord((0, 0))
            >>> Chord((0, 0, 1, 1, 2, 2)).remove_chord(0)
            Chord((0, 0, 1, 1))
            >>> Chord((0, 1, 1, 0, 2, 3, 3, 2)).remove_chord(2)
            Chord((0, 1, 1, 0, 2, 2))
        """
        if selected is None:
            if self._length == 0:
                return self
            selected = self._length - 1
        assert 0 <= selected and selected < self._length
        return Chord(
            list(val if val < selected else val - 1 for val in self if val != selected)
        )

    # sCP, sAsk: it feels like there is a better way to do this
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

    # sToDo: more tests. sCP: similar method for this was implemented in occurances_in
    def _contains(self, patt: "Chord") -> bool:
        """Determines if a chord contains an instance of patt"""
        if isinstance(patt, Chord):
            return len(patt.occurrences_in(self)) > 0
        raise TypeError("patt must be a Chord")
    
    # sCP, sToDo: more tests, changed function in structure (just calls occurrences_in)
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
                subchord_indices_in_patt.extend(patt.chord[i])
            subchord_indices_in_patt.sort()

            subchord_patt = []
            for i in subchord_indices_in_patt:
                subchord_patt.append(patt.get_chord()[i])

            if Chord(()).standardize(subchord_patt) == self:
                append = True
                if patt_has_colour and self_has_colour:
                    for i in range(len(subchord_indices_in_patt)):
                        if colour_patt[subchord_indices_in_patt[i]] != colour_self[i]:
                            append = False
                if append:
                    instances.append(subchord)

        return instances

    # not really needed with _pattern parameter. OR make default in get_chords method
    def get_chord(self) -> "Chord":
        """Returns the chord part of the pattern.

        Examples:
            >>> Chord((3,2,1,0)).get_chord()
            Chord((3, 2, 1, 0))
        """
        return self

    # sCP
    def avoids(self, *chords: "Chord") -> bool:
        """
        Check if self avoids all chords.
        """
        return all(not self._contains(patt) for patt in chords)

    # sCP, some changes to arguement format: m,n -> (m,n)
    def insert(self, pos: Optional[Iterable[int]] = None, m: Optional[int] = None, n: Optional[int] = None) -> "Chord":
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
        if pos == None or len(pos) !=2:
            return self.concat(Chord((0, 0)))
        
        m = pos[0]
        n = pos[1]
        
        assert (0 <= m) and (m <= n) and (n <= listlen)
        result = list(self)
        result.insert(m, -1)
        result.insert(n+1, -1)
        return self._standardize(result)

    # sCP    
    def get_chords(self, selected: list) -> "Chord":
        """Return the Chord based on the selected list of foots.

        Examples:
            >>> Chord((0, 0, 1, 1)).get_chords([0])
            Chord((0, 0))
            >>> Chord((0, 1, 0, 1)).get_chords([0, 1])
            Chord((0, 1, 0, 1))
        """
        return self.standardize([val for val in self if val in selected])

    # sCP, changed name from chord -> chord_dict
    @property
    def chord(self) -> dict[(int)]:
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

class GriddedChord(CombinatorialObject):
    def __init__(
        self, pattern: Iterable[int] = (), positions: Iterable[Cell] = ()
    ) -> None:
        self._patt = Chord(pattern)
        self._pos = tuple(positions)
        if len(self._patt) * 2 != len(self._pos):
            raise ValueError("Pattern and position list have unequal lengths.")
        self._cells: FrozenSet[Cell] = frozenset(self._pos)

    @classmethod
    def single_cell(cls, pattern: Iterable[int], cell: Cell) -> "GriddedChord":
        """Construct a gridded chord where the chords are all located in a
        single cell."""
        return cls(pattern, (cell for _ in range(int(len(pattern)))))

    @classmethod
    def empty_chord(cls) -> "GriddedChord":
        """Construct the empty gridded chord."""
        return cls((), ())

    @classmethod 
    # sCN: point_chord -> single_chord, sDone
    def single_chord(cls, cells: Iterable[Cell]) -> "GriddedChord":
        """Construct the single gridded chord using the cells given. If only one cell 
        is given, will construct the single gridded chord with both ends in the same cell"""
        if len(cells) == 1:
            return cls((0,0), (cells[0],cells[0]))
        elif len(cells) == 2:
            return cls((0,0), (cells[0],cells[1]))
        raise ValueError("incorrect number of cells given")

    def occupies(self, cell: Cell) -> bool:
        """Checks if the gridded chord has a point in the given cell."""
        return cell in self._cells

    # sToDo: could be generalized to recognize patterns that are not indexed the same
    def occurrences_in(self, other: "GriddedChord") -> Iterator[Tuple[int, ...]]:
        """Returns all occurrences of self in other."""
        yield from self._patt.occurrences_in(other._patt, self._pos, other._pos)

    def occurs_in(self, other: "GriddedChord") -> bool:
        """Checks if self occurs in other."""
        return any(self.occurrences_in(other))

    def avoids(self, *patts: "GriddedChord") -> bool:
        """Return true if self avoids all of the patts."""
        return not self.contains(*patts)

    def contains(self, *patts: "GriddedChord") -> bool:
        """Return true if self contains an occurrence of any of patts."""
        return any(any(True for _ in patt.occurrences_in(self)) for patt in patts)
    
    #sCN: the old method was perm specific
    def contradictory(self) -> bool:
        """Checks if the points of the griddedchord contradict the chord.

        Checks if for every 0 <= i < j < n, 
            if i and j are on the same chord, their x positions are the same
            if one chord comes before the other, the x position is less or equal
            the y positions are in increasing order 
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
        chord."""
        cells = set(cells)
        chords_to_remove = [chord for chord, pos in self if pos in cells]
        remaining_chords = []
        remaining_positions = []

        for i, (chord, pos) in enumerate(self):
            if chord not in chords_to_remove:
                remaining_chords.append(chord)
                remaining_positions.append(pos)

        return type(self)(
            Chord(()).standardize(remaining_chords)._pattern,
            remaining_positions
        )

    def points_in_cell(self, cell: Cell) -> Iterator[int]:
        """Yields the indices of the points in the cell given."""
        return (i for i, (_, pos) in enumerate(self) if pos == cell)

    def isolated_cells(self) -> Iterator[Cell]:
        """Yields the cells that contain only one chord of the gridded
        chord and are in their own row and column (besides the other end of the chord)."""
        return (
            (x1, y1)
            for i, (chord1, (x1, y1)) in enumerate(self)
            if not any(
                x1 == x2 or y1 == y2 for j, (chord2, (x2, y2)) in enumerate(self) 
                if (i != j and chord1 != chord2)
            )
        )

    # sCN: added check for chords being in same column
    def is_isolated(self, indices: Iterable[int]) -> bool:
        """Checks if the cells at the indices do not share a row or column with
        any other cell in the gridded chord."""
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
        respect to the given force."""
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
        """Yields all points of the gridded chord in the column col."""
        return ((val, i) for i, (val, (x, _)) in enumerate(self) if x == col)

    def get_points_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points of the gridded chord in the row."""
        return ((val, i) for i, (val, (_, y)) in enumerate(self) if y == row)

    def get_points_below_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points of the gridded chord below the row."""
        return ((val, i) for i, (val, (_, y)) in enumerate(self) if y < row)

    def get_points_above_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points of the gridded chord above the row."""
        return ((val, i) for i, (val, (_, y)) in enumerate(self) if y > row)

    def get_points_left_col(self, col) -> Iterator[Tuple[int, int]]:
        """Yields all points of the gridded chord left of column col."""
        return ((val, i) for i, (val, (x, _)) in enumerate(self) if x < col)
    
    def get_points_right_col(self, col: int) -> Iterator[Tuple[int, int]]:
        """Yields all points of the gridded chord right of column col."""
        return ((val, i) for i, (val, (x, _)) in enumerate(self) if x > col)

    def get_subchord_left_col(self, col: int) -> "GriddedChord":
        """Returns the gridded subchord of points left of column col."""
        gen1, gen2 = tee((i for chord, i in self.get_points_left_col(col)), 2)
        return type(self)(
            Chord(()).standardize([self._patt[i] for i in gen1])._pattern, (self._pos[i] for i in gen2)
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
            Chord(()).standardize([self._patt[i] for i in gen1])._pattern,
            (self._pos[i] for i in gen2),
        )

    # sCN: used cleaner list notation
    def remove_chord(self, index) -> "GriddedChord":
        avoid_chord = self._patt[index]
        new_chords = [chord for chord, _ in self if chord != avoid_chord]
        new_positions = [pos for chord, pos in self if chord != avoid_chord]

        new_patt = Chord(()).standardize(new_chords)._pattern

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

        new_grid = GriddedChord(Chord(()).standardize(patt)._pattern, positions)


        return new_grid

    # sCN: previous method returned the restrictions on the values for a point inserted 
    # into a given cell for permutations. Changed to accomadate chords.
    def get_bounding_box(self, col: int, row1: int, row2: int) -> Tuple[int, int, int, int]:
        """Determines the range possible indices for a chord inserted into column col, with bottom 
        endpoint in row1, and top endpoint in row2. Returns in order: (min chord, max chord)"""
        row = list(self.get_points_row(cell[1]))
        cols_left = list(self.get_points_col(cell[0]))
        # row is empty
        if not row:
            above = list(self.get_points_above_row(cell[1]))
            # above is empty
            if not above:
                maxval = len(self)
            else:
                maxval = min(p[1] for p in above)
            below = list(self.get_points_below_row(cell[1]))
            if not below:
                minval = 0
            else:
                minval = max(p[1] for p in below) + 1
        else:
            maxval = max(p[1] for p in row) + 1
            minval = min(p[1] for p in row)
        if not col:
            right = list(self.get_points_right_col(cell[0]))
            if not right:
                maxdex = len(self)
            else:
                maxdex = min(p[0] for p in right)
            left = list(self.get_points_left_col(cell[0]))
            if not left:
                mindex = 0
            else:
                mindex = max(p[0] for p in left) + 1
        else:
            maxdex = max(p[0] for p in col) + 1
            mindex = min(p[0] for p in col)
        return (minchord, maxchord, mindex, maxdex)
    '''
    def insert_point(self, cell: Cell) -> Iterator["GriddedChord"]:
        """Insert a new point into cell of the gridded chord, such that the
        point is added to the underlying pattern with the position at the cell.
        Yields all gridded chords where the point has been mixed into the points
        in the cell."""
        mindex, maxdex, minval, maxval = self.get_bounding_box(cell)
        for idx, val in product(range(mindex, maxdex + 1), range(minval, maxval + 1)):
            yield self.insert_specific_point(cell, idx, val)

    def insert_specific_point(self, cell: Cell, idx: int, val: int) -> "GriddedChord":
        """Insert a point in the given cell with the given idx and val."""
        return type(self)(
            self._patt.insert(idx, val), self._pos[:idx] + (cell,) + self._pos[idx:]
        )

    def remove_point(self, index: int) -> "GriddedChord":
        """Remove the point at index from the gridded chord."""
        patt = Chord.to_standard(self.patt[:index] + self.patt[index + 1 :])
        pos = self.pos[:index] + self.pos[index + 1 :]
        return type(self)(patt, pos)

    def all_subchords(self, proper: bool = True) -> Iterator["GriddedChord"]:
        """Yields all gridded subchords."""
        for r in range(len(self) if proper else len(self) + 1):
            for subidx in combinations(range(len(self)), r):
                yield type(self)(
                    Chord.to_standard(self._patt[i] for i in subidx),
                    (self._pos[i] for i in subidx),
                )

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
                self.insert_specific_point((x, y), index, n)
                for index in range(n + 1)
                for x in range(
                    self.pos[index - 1][0] if index > 0 else 0,
                    self.pos[index][0] + 1 if index < n else c,
                )
                for y in range(min_y, r)
            )

    def apply_map(self, cell_mapping: Callable[[Cell], Cell]) -> "GriddedChord":
        """Map the coordinates to a new list of coordinates according to the
        cell_mapping given."""
        return type(self)(self._patt, [cell_mapping(cell) for cell in self._pos])

    def is_point_chord(self) -> bool:
        """Checks if the gridded chord is of length 1."""
        return len(self) == 1

    def is_localized(self) -> bool:
        """Check if the gridded chord occupies only a single cell."""
        return self.is_single_cell()

    def is_single_cell(self) -> bool:
        """Check if the gridded chord occupies only a single cell."""
        return len(self._cells) == 1

    def is_single_row(self) -> bool:
        """Check if the gridded chord occupies only a single row."""
        return len(set(y for (_, y) in self._cells)) == 1

    def is_empty(self) -> bool:
        """Check if the gridded chord is the gridded chord."""
        return not bool(self._patt)

    def is_interleaving(self) -> bool:
        """Check if the gridded chord occupies two cells that are in the
        same row or column."""
        xs: Dict[int, int] = {}
        ys: Dict[int, int] = {}
        for x, y in self._cells:
            if xs.get(x, y) != y or ys.get(y, x) != x:
                return True
            xs[x] = y
            ys[y] = x
        return False

    def factors(self) -> List["GriddedChord"]:
        """Return a list containing the factors of a gridded chord.
        A factor is a sub gridded chord that is isolated on its own rows
        and columns."""
        uf = UnionFind(len(self.pos))
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
    def decompress(cls, arr: array) -> "GriddedChord":
        """Decompresses a list of integers in the form outputted by the
        compress method and constructs an Obstruction."""
        n, it = len(arr) // 3, iter(arr)
        return cls(
            Chord(next(it) for _ in range(n)), ((next(it), next(it)) for _ in range(n))
        )

    # Symmetries
    def reverse(self, transf: Callable[[Cell], Cell]) -> "GriddedChord":
        """
        Reverses the tiling within its boundary. Every cell and obstruction
        gets flipped over the vertical middle axis."""
        return type(self)(self._patt.reverse(), reversed(list(map(transf, self._pos))))

    def complement(self, transf: Callable[[Cell], Cell]) -> "GriddedChord":
        """Flip over the horizontal axis."""
        return type(self)(self._patt.complement(), map(transf, self._pos))

    def inverse(self, transf: Callable[[Cell], Cell]) -> "GriddedChord":
        """Flip over the diagonal."""
        flipped = self._patt.inverse()
        pos = self._patt.inverse().apply(self._pos)
        return type(self)(flipped, map(transf, pos))

    def antidiagonal(self, transf: Callable[[Cell], Cell]) -> "GriddedChord":
        """Flip over the anti-diagonal"""
        flipped = self._patt.flip_antidiagonal()
        pos = self._patt.rotate(-1).apply(self._pos)
        return type(self)(flipped, map(transf, pos))

    def column_reverse(self, column: int) -> "GriddedChord":
        """Reverse the part of the gridded chord that belongs to a given column."""
        parts: Tuple[List[Tuple[int, Tuple[int, int]]], ...] = ([], [], [])
        for v, (x, y) in self:
            parts[(x >= column) + (x > column)].append((v, (x, y)))
        return type(self)(*zip(*chain(parts[0], reversed(parts[1]), parts[2])))

    def row_complement(self, row: int) -> "GriddedChord":
        """Replace the part of the gridded chord that belongs to a row with its
        complement."""
        indices, vals = set(), []
        for i, (val, (_, y)) in enumerate(self):
            if y == row:
                indices.add(i)
                vals.append(val)
        st = Chord.to_standard(vals)
        unstandardized = dict(zip(st, vals))
        in_row = (unstandardized[val] for val in st.complement())
        return type(self)(
            Chord(
                next(in_row) if i in indices else val for i, val in enumerate(self.patt)
            ),
            self.pos,
        )

    def chordute_columns(self, chord: Iterable[int]) -> "GriddedChord":
        """Given an initial state of columns 12...n, chordute them using the provided
        chord.
        """
        if not isinstance(chord, Chord):
            chord = Chord(chord)
        assert len(chord) > max(x for x, y in self.pos)
        cols: List[List[Tuple[int, int]]] = [[] for _ in range(len(chord))]
        for v, (x, y) in self:
            cols[x].append((v, y))
        patt, positions = [], []
        for x, col in enumerate(chord.apply(cols)):
            for v, y in col:
                patt.append(v)
                positions.append((x, y))
        return type(self)(patt, positions)

    def chordute_rows(self, chord: Iterable[int]) -> "GriddedChord":
        """Given an initial state of rows 12...n, chordute them using the provided
        chord.
        """
        if not isinstance(chord, Chord):
            chord = Chord(chord)
        assert len(chord) > max(x for x, y in self.pos)
        back_map = dict(zip(chord, range(len(chord))))
        positions = [(x, back_map[y]) for x, y in self.pos]
        occ: List[List[int]] = [[] for _ in range(len(chord))]
        for val, (_, y) in zip(self.patt, positions):
            occ[y].append(val)
        offset, val_map = 0, {}
        for lis in occ:
            for a, b in zip(lis, Chord.to_standard(lis)):
                val_map[a] = b + offset
            offset += len(lis)
        return type(self)((val_map[val] for val in self.patt), positions)

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
        max_x = max(cell[0] for cell in self.pos)
        max_y = max(cell[1] for cell in self.pos)
        res = ""

        def points_in_col(i):
            return sum(1 for cell in self.pos if cell[0] == i)

        def points_in_row(j):
            return sum(1 for cell in self.pos if cell[1] == j)

        row_boundary = (
            "+" + "+".join("-" * points_in_col(i) for i in range(max_x + 1)) + "+"
        )
        col_boundary = (
            "|" + "|".join(" " * points_in_col(i) for i in range(max_x + 1)) + "|"
        )

        for j in range(max_y, -1, -1):
            res += (
                "\n".join(
                    [row_boundary] + [col_boundary for i in range(points_in_row(j))]
                )
                + "\n"
            )
        res += row_boundary

        for (idx, val) in enumerate(self.patt):
            x, y = self.pos[idx]
            # insert into this spot:
            # (idx + x + 1) is the horizontal index. idx is points to left, and
            #               x + 1 counts number of col boundaries to left
            # (len(self) + max_y) - (val + y) is vertical index.
            #               val is points below, and y + 1 counts number of - below
            insert = (idx + x + 1) + ((len(self) + max_y) - (val + y)) * (
                len(col_boundary) + 1
            )
            res = res[:insert] + "●" + res[insert + 1 :]
        return res
'''

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

gc = GriddedChord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4), ((0, 0), (1, 0), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4)))
hc = GriddedChord((0, 1, 0, 2,2,1), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1) ))
print(gc.get_gridded_chord_in_cells(((0,0), (1,1))))
