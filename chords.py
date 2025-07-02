from comb_spec_searcher import CombinatorialObject
from typing import Callable, Dict, FrozenSet, Iterable, Iterator, List, Optional, Tuple, Union
from itertools import combinations, islice, tee, product, chain, filterfalse
from permuta.misc import UnionFind
import json
from collections import deque

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
    def __init__(self, iterable: Iterable[int] = ()) -> None:
        # Cache for data used when finding occurrences of self in a perm
        self._cached_pattern_details: Optional[List[Tuple[int, int, int, int]]] = None
        self._length = len(iterable)//2
        self._pattern = tuple(iterable)
        self._patt_length = len(iterable)

        # Creates dictionary of chord number and vertices it connects
        # (0, 1, 2, 1, 3, 3, 2, 0) becomes {0: (0, 7), 1: (1, 3), 2: (2, 6), (3: (4, 5)}
        self._chord_dict = {}
        for i,c in enumerate(iterable):
            if c not in self._chord_dict:
                self._chord_dict[c] = (i, -1)
            else:
                self._chord_dict[c] = (self._chord_dict[c][0],i)

    def assert_valid_chord(self) -> None:
        """Examples:
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
        TypeError: 'None' object is not an integer"""
        iterable = self._pattern

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
            return []
        
        if self_has_colour and len(colour_self) != self._patt_length:
            raise Exception("Incorrect colour length for self")
        
        if patt_has_colour and len(colour_patt) != patt._patt_length:
            raise Exception("Incorrect colour length for pattern given")

        if len(self) > len(patt):
            return []
        
        instances = []
        
        indexed_patt = list(enumerate(patt._pattern)) # matches each val in patt to an index
        self_reindexed = Chord.reindex(self._pattern) # makes sure self is in "normal" form

        # loops over patterns of length self in patt
        for subslice in combinations(indexed_patt, self._patt_length): 
            subpatt = [val for _, val in subslice] # extracts pattern information from sub pattern
            if Chord.reindex(subpatt) == self_reindexed:
                append = True
                patt_indices = [idx for idx, _ in subslice] # extractes index information about sub pattern
                if patt_has_colour and self_has_colour:
                    # loops over indices of self and subpattern over length of patterns
                    for idx_self, idx_patt in enumerate(patt_indices): 
                        if colour_patt[idx_patt] != colour_self[idx_self]: 
                            append = False 
                if append:
                    instances.append(tuple(set(subpatt)))

        return instances

    @classmethod
    def reindex(cls, patt: Iterable) -> "Chord":
        patt_vals = list(set(patt))
        patt_vals.sort()
        chord_map = {}

        for i, val in enumerate(patt_vals):
            chord_map[val] = i

        reindexed_patt = []
        for i in patt:
            reindexed_patt.append(chord_map[i])

        return Chord(reindexed_patt)

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
        return self._chord_dict
    
    # sToDo: test
    @classmethod
    def __generate(cls, n: int):
        """
        generate all chord diagrams with size n
        """
        length = 2*n
        q = deque()
        q.append(((1 << length) - 1, [-1]*length, 0))
        while q:
            remaining, chord, i = q.popleft()
            if remaining == 0:
                yield cls(chord)
                continue
            source = 0
            while not (remaining & (1 << (length - source-1))):
                source += 1
            remaining ^= (1 << (length - source-1))
            for sink in range(source+1, length):
                if remaining & (1 << (length - sink - 1)):
                    new_chord = chord[:]
                    new_chord[source] = i
                    new_chord[sink] = i
                    new_remaining = remaining^(1 << (length - sink - 1))
                    q.append((new_remaining, new_chord, i+1))

    # sToDo: test
    @classmethod
    def of_length(cls, length: int) -> Iterator["Chord"]:
        """Generate all chord diagrams of a given length.

        Examples:
            >>> list(Chord.of_length(2))
            [Chord(0,1,1,0),Chord((0, 1, 0, 1)), Chord((0, 0, 1, 1))]
        """
        yield from cls.__generate(length)

    # sToDo: all of the following up to connected (these methods will probably be fun)
    def crossing(self):
        pass
    def nesting(self):
        pass
    def k_crossing(self):
        pass

    def connected(self):
        nodes = list(range(len(self)))
        edges = Chord((0, 1, 0, 1)).occurrences_in(self)
        graph = [[] for _ in nodes]
        for edge in edges:
            u = edge[0]
            v = edge[1]
            graph[u].append(v)
            graph[v].append(u)

        visited = [False for _ in nodes]

        # Depth first search for graph to see what is connected
        def dfs(node, graph, visited):
            visited[node] = True
            # search all unvisited neighbors
            for neighbor in graph[node]:
                if not visited[neighbor]:
                    dfs(neighbor, graph, visited)
        if len(graph) > 0:
            dfs(0, graph, visited)
        else: # empty graph is disconnected
            return False

        #print(nodes, edges, graph, visited, all(visited))

        return all(visited)

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

        # Creates dictionary of chord number and vertices it connects
        # (0, 1, 2, 1, 3, 3, 2, 0) becomes {0: (0, 7), 1: (1, 3), 2: (2, 6), (3: (4, 5)}
        self._chord_dict = self._chord.chord_dict()

        # this was needed, but has problem with point placement when you try to place a single point and then another one
        #if self.contradictory():
        #    raise ValueError("Contradictory positions given.")
    
    ### Generators for new GriddedChords ###
    @classmethod
    def empty_chord(cls) -> "GriddedChord":
        """Construct the empty gridded chord.

        Examples:
        GriddedChord.empty_chord() == GriddedChord()
        GriddedChord.empty_chord() == GriddedChord(Chord(), ())"""
        return cls(Chord(), ())
    
    @classmethod 
    def single_chord(cls, cells: Iterable[Cell]) -> "GriddedChord":
        # sCN: point_chord -> single_chord
        """Construct the single gridded chord using the cells given. If only one cell 
        is given, will construct the single gridded chord with both ends in the same cell
        
        Examples:
        >>> GriddedChord.single_chord((0,0))
        GriddedChord((0, 0), ((0, 0), (0, 0)))
        >>> GriddedChord.single_chord((0, 1), (1, 1))
        GriddedChord((0, 0), ((0, 1), (1, 1)))"""
        if len(cells) == 1:
            return cls(Chord((0,0)), (cells[0],cells[0]))
        elif len(cells) == 2:
            if cells[0][1] != cells[1][1]:
                raise Exception("inconsistant rows given")
            return cls(Chord((0,0)), (cells[0],cells[1]))
        raise ValueError("incorrect number of cells given")
   
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
    def all_grids(cls, chord: Chord, cells) -> Iterator["GriddedChord"]:
        """Returns all possible griddings of chord in cells
        
        Examples: 
        >>> GriddedChord.all_grids(Chord((0, 0, 1, 1)), [(0, 0)])
        GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))
        >>> GriddedChord.all_grids(Chord((0, 0, 1, 1)), [(0, 0), (0, 1)])
        GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 1), (0, 1)))
        GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 1), (0, 1)))
        GriddedChord(Chord((0, 0, 1, 1)), ((0, 0), (0, 0), (0, 0), (0, 0)))
        """
        positions_so_far = [[]]
        source_vals = {} # chord value : source
        last_source_idx = -1
        
        def append_pos_lists(curr_positions: List[List[Cell]], 
                        possible_cells: List[Cell], 
                        source_idx: int = None, 
                        last_source_idx: int= None):
            """Returns the list of all positions that can be created by extending a
            position in positions_curr with a cell in positions_possible
            
            source_idx is the source index of the chord in the position about to be added
            
            last_source_idx is the index of where the largest chord so far started"""
            appended_positions = []
            # loops through every current position list, and checks (then adds) every possible cell as a new list
            for pos_lst in curr_positions:
                for cell in possible_cells:
                    # x-coordinate of the most recent cell
                    last_x = pos_lst[-1][0] if len(pos_lst) > 0 else -1
                    valid_x = (cell[0] >= last_x) # the x-coordinates must be in increasing order

                    valid_y = True
                    # if the chord has a source, the only valid y-coordinate is the same as the source
                    if source_idx != None:
                        valid_y = (pos_lst[source_idx][1] == cell[1]) 
                    # if a new (non-root) chord is being placed, it must go above the last chord placed
                    elif (last_source_idx >= 0):
                        valid_y = (pos_lst[last_source_idx][1] <= cell[1])
                    
                    # extend pos_lst by cell iff it has a valid x and y coordinate.
                    if valid_x and valid_y:
                        pos_lst_extended = pos_lst.copy()
                        pos_lst_extended.append(cell)
                        appended_positions.append(pos_lst_extended)
            return appended_positions

        # runs over every item in the chord to build the posible position lists to the length of the chord
        for idx, chord_val in enumerate(chord):
            source = source_vals.get(chord_val) # either the idx of the source or none if chord_val is a source
            # updates the possible positions lists to include positions for (idx, chord_val)
            positions_so_far = append_pos_lists(positions_so_far, cells, source, last_source_idx)
            if source == None: # if source is none, the current chord is the most recent source
                source_vals[chord_val] = idx
                last_source_idx = idx
            
        for pos_list in positions_so_far:
            if (len(pos_list) == len(chord) * 2):
                yield cls(chord, pos_list)

    @classmethod
    def point_chord(cls, pos: Cell) -> "GriddedChord":
        return GriddedChord(Chord((0,)), pos)

    ### Return information about GriddedChord instance ###
    def occurrences_in(self, other: "GriddedChord") -> Iterator[Tuple[int, ...]]:
        # sToDo: could be generalized to recognize patterns that are not indexed the same
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

    def points_in_cell(self, cell: Cell) -> Iterator[int]:
        """Yields the indices of the chords in the cell given.
        
        Examples: 
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1))).points_in_cell((1, 1))
        (0, 1, 2, 3, 4, 5)
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0, 1), (1, 1), (0, 1), (1, 1), (1, 1), (1, 1))).points_in_cell((1, 1))
        (2, 3, 4, 5)"""
        return (i for i, (_, pos) in enumerate(self) if pos == cell)
    
    def in_cell(self, cell: Cell) -> bool:
        return any(pos == cell for i, (_, pos) in enumerate(self))

    def isolated_cells(self) -> Iterator[Cell]:
        """Yields the cells that contain only one chord of the gridded
        chord and are in their own row and column (besides the other end of the chord).
        """
        return set(
            (x1, y1)
            for i, (chord1, (x1, y1)) in enumerate(self)
            if not any(
                x1 == x2 or y1 == y2 for j, (chord2, (x2, y2)) in enumerate(self) 
                if (i != j and chord1 != chord2))
        )

    def forced_point_index(self, cell: Cell, direction: int) -> int:
    # sCN: returns most extreme chord in each case, not point
        """Search in the cell given for the chord with the strongest force with
        respect to the given force.
        
        Examples:
        >>> GriddedChord(Chord((0, 1, 0, 1, 2, 2)), ((0,0), (0,0), (0,0), (0,0), (0, 0), (0, 0))).forced_point_index((0, 0), DIR_EAST)
        2
        >>> GriddedChord(Chord((0, 1, 0, 2, 2, 1)), ((0,0), (0,0), (0,0), (0,0), (0, 0), (0, 0))).forced_point_index((0, 0), DIR_NORTH)
        1
        >>> GriddedChord(Chord((0, 1, 0, 1, 2, 2)), ((0,0), (0,0), (0,0), (0,0), (0, 0), (0, 0))).forced_point_index((0, 0), DIR_WEST)
        0
        >>> GriddedChord(Chord((0, 1, 0, 1, 2, 2)), ((0,0), (0,0), (0,0), (0,0), (0, 0), (0, 0))).forced_point_index((0, 0), DIR_SOUTH)
        0"""
        if self.occupies(cell):
            indices = self.points_in_cell(cell)
            chords = (chord for chord, cell in self if cell == cell)
            if direction == DIR_EAST:
                return self._patt[max(indices)]
            if direction == DIR_NORTH:
                return max(chords)
            if direction == DIR_WEST:
                return self._patt[min(indices)]
            if direction == DIR_SOUTH:
                return min(chords)
            raise ValueError("You're lost, no valid direction")
        raise ValueError("The gridded chord does not occupy the cell")

    def chord_dict(self) -> Dict[(int)]:
        return self._chord_dict

    ### Retrun Boolean information ###
    def is_empty(self) -> bool:
        """Check if the gridded chord is the empty gridded chord."""
        return len(self) == 0
    
    def is_point(self) -> bool:
        return self._chord._patt_length == 1

    def is_single_chord(self) -> bool:
        # sCN: is_point_chord -> is_single_chord
        """Checks if the gridded chord is of length 1."""
        return len(self) == 1

    def is_localized(self) -> bool:
        # sToDo: why is this function just calling another? redundant?
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

    def is_isolated(self, indices: Iterable[int]) -> bool:
        # sCN: added check for chords being in same column
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

    def is_interleaving(self) -> bool:
        # sToDo: change to ignore chords in same column?
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
    
    def is_connected(self, cells: List[Cell] = []) -> bool:
        """Determines whether a diagram is connected.
        
        cells: list of cells to check if the subchord inside is connected.
                if none are passed, the whole chord will be checked.
                
        returns: whether the chord based on cells is connected
        
        Examples:
        >>> GriddedChord(Chord((0, 1, 1, 0), ((0, 0), (1, 0), (1, 0), (2, 0))).is_connected())
        False
        >>> GriddedChord(Chord((0, 1, 0, 1), ((0, 0), (1, 0), (0, 0), (2, 0))).is_connected())
        True
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 0), (2, 0))).is_connected([(0, 0)]))
        True"""
        # sAsk: is linkage everything strictly inside, or everything with one end inside. (strictly for now)
        subchord = self
        if len(cells) != 0:
            subchord = self.get_subchord_in_cells(cells)

        return subchord._chord.connected()
 
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
        >>> gc = GriddedChord(Chord((0, 1, 2, 0, 3, 4, 2, 3, 1, 4)), 
                               ((0, 0), (1, 0), (1, 1), (0, 1), (3, 2), (4, 2), (1, 2), (3, 3), (1, 3), (4, 4)))
        >>> gc.contains(GriddedChord(Chord((0, 1, 0, 1)), ((0, 0), (1, 0), (0, 1), (1, 3))))
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
                if (chordj == chordi) and (yj != yi):
                    return True
                
                # The x positions of increasing chords must increase (or stay the same)
                if ((chordj < chordi) and (yj > yi)) or ((chordj > chordi) and (yj < yi)):
                    return True
                    
                # y positions must be in increasing order of the chord
                if (xj < xi):
                    return True
                    
        return False
        

    ### Get set of points or subchord ###
    def get_points_col(self, col: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord in the column col.
        
        Examples:
        >>> gc=GriddedChord(Chord((0,1,2,0,3,4,2,3,1,4)), ((0,0), (1,0), (1,1), (0,1), (3,2), (4,2), (1,2), (3,3), (1,3), (4,4)))
        >>> gc.get_points_col(1)
        (1, 1), (1, 8), (2, 2), (2, 6)
        >>> gc.get_points_col(2)
        ()
        >>> gc.get_points_col(3)
        (3, 4), (3, 7)"""
        return ((i, val) for i, (val, (x, _)) in enumerate(self) if x == col)

    def get_points_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord in the row.
        
        Examples:
        >>> gc=GriddedChord(Chord((0,1,2,0,3,4,2,3,1,4)), ((0,0), (1,0), (1,1), (0,1), (3,2), (4,2), (1,2), (3,3), (1,3), (4,4)))
        >>> gc.get_points_row(0)
        (0, 0), (1, 1)
        >>> gc.get_points_row(5)
        ()
        >>> gc.get_points_row(1)
        (0, 3), (2, 2)"""
        return ((i, val) for i, (val, (_, y)) in enumerate(self) if y == row)

    def get_points_below_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord below the row.
        
        Examples:
        >>> gc=GriddedChord(Chord((0,1,2,0,3,4,2,3,1,4)), ((0,0), (1,0), (1,1), (0,1), (3,2), (4,2), (1,2), (3,3), (1,3), (4,4)))
        >>> gc.get_points_below_row(2)
        (0, 0), (0, 3), (1, 1), (2, 2)
        >>> gc.get_points_below_row(0)
        ()
        """
        return ((i, val) for i, (val, (_, y)) in enumerate(self) if y < row)

    def get_points_above_row(self, row: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord above the row.
        
        Examples:
        >>> gc=GriddedChord(Chord((0,1,2,0,3,4,2,3,1,4)), ((0,0), (1,0), (1,1), (0,1), (3,2), (4,2), (1,2), (3,3), (1,3), (4,4)))
        >>> gc.get_points_above_row(2)
        (1, 8), (3, 7), (4, 9)
        >>> gc.get_points_above_row(4)
        ()
        >>> 
        """
        return ((i, val) for i, (val, (_, y)) in enumerate(self) if y > row)

    def get_points_left_col(self, col) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord left of column col.
        
        Examples:
        >>> gc=GriddedChord(Chord((0,1,2,0,3,4,2,3,1,4)), ((0,0), (1,0), (1,1), (0,1), (3,2), (4,2), (1,2), (3,3), (1,3), (4,4)))
        >>> gc.get_points_left_col(0)
        ()
        >>> gc.get_points_left_col(1)
        (0, 0), (0, 3)
        >>> gc.get_points_left_col(2)
        (0, 0), (0,3 ), (1, 1), (1, 8), (2, 2), (2, 6)
        """
        return ((i, val) for i, (val, (x, _)) in enumerate(self) if x < col)
    
    def get_points_right_col(self, col: int) -> Iterator[Tuple[int, int]]:
        """Yields all points in form (chord, index) of the gridded chord right of column col.
        
        Examples
        >>> gc=GriddedChord(Chord((0,1,2,0,3,4,2,3,1,4)), ((0,0), (1,0), (1,1), (0,1), (3,2), (4,2), (1,2), (3,3), (1,3), (4,4)))
        >>> gc.get_points_right_col(4)
        ()
        >>> gc.get_points_right_col(3)
        (4, 5), (4, 9)
        >>> gc.get_points_right_col(2)
        (4, 5), (4, 9), (3, 4), (3, 7)
        """
        return ((i, val) for i, (val, (x, _)) in enumerate(self) if x > col)

    def get_subchord_below_row(self, row: int) -> "GriddedChord":
        """Returns the gridded subchord of points left of column col.
        
        Examples:
        gc=GriddedChord(Chord((0,1,2,0,3,4,2,3,1,4)), ((0,0), (1,0), (1,1), (0,1), (3,2), (4,2), (1,2), (3,3), (1,3), (4,4)))
        >>> gc.get_subchord_left_col(0)
        GriddedChord()
        >>> gc.get_subchord_left_col(1)
        GriddedChord(Chord(0, 0), ((0, 0), (0, 1)))
        >>> gc.get_subchord_left_col(2)
        GriddedChord(Chord(0, 1, 2, 0, 2, 1), ((0, 0), (1, 0), (1, 1), (0, 1), (1, 2), (1, 3)))"""
        return type(self)(
            Chord.to_standard([chord for _, chord in self.get_points_below_row(row)]), 
            (self._pos[i] for i, _ in self.get_points_below_row(row))
        )
    
    def get_subgrid_at_chords(self, chords: Iterable[int]) -> "GriddedChord":
        # sCN: used to be from indices, now from chords 
        # (previously called get_gridded_chord_at_indices(self, indices: Iterable[int]))
        """ Returns the subgridded chord that contains only the given chords"""
        indices = [idx for idx, (chord, pos) in enumerate(self) if chord in chords]
        return type(self)(
            Chord.to_standard([self._patt[i] for i in indices]),
            (self._pos[i] for i in indices),
        )
    
    def get_subchord_in_cells(self, cells: Iterable[Cell]) -> "GriddedChord":
        """Returns the subgridded chord with chords with both endpoints in cells.
        """
        chords_to_throw = [chord for chord, pos in self if pos not in cells]
        patt = []
        positions = []

        # builds lists of the chord pattern and positions of the chords in cells
        for chord, pos in self:
            if chord not in chords_to_throw:
                patt.append(chord)
                positions.append(pos)

        new_grid = GriddedChord(Chord.to_standard(patt), positions)
        return new_grid
    
    def get_chords_on_cells(self, cells: Iterable[Cell]) -> "GriddedChord":
        """Returns the subgridded chord with chords with an endpoint in cells.
        """
        chords_to_keep = [chord for chord, pos in self if pos in cells]
        patt = []
        positions = []

        # builds lists of the chord pattern and positions of the chords in cells
        for chord, pos in self:
            if chord in chords_to_keep:
                patt.append(chord)
                positions.append(pos)

        new_grid = GriddedChord(Chord.to_standard(patt), positions)
        return new_grid

    def all_subchords(self, proper: bool = True) -> Iterator["GriddedChord"]:
        """Yields all gridded subchords."""
        # loops through all sizes of subchords, taking into account if we are looking for proper subchords or not
        for subchord_length in range(len(self) if proper else len(self) + 1):
            for subchords in combinations(range(len(self)), subchord_length):
                yield type(self)(
                    Chord.to_standard([chord for chord, pos in self if chord in subchords]),
                    (pos for chord, pos in self if chord in subchords))

    def get_bounding_indices(self, col_source: int, col_sink: int) -> Tuple[int, int, int, int]:
        # sCN: previous method returned the restrictions on the values for a point inserted 
        # into a given cell for permutations. Changed to accomadate chords.
        """Determines the range possible indices where a chord with source in col_source and sink 
        in col_sink can be inserted into self.
        Returns in order: (min index source, max index source, min index sink, max index sink)
        
        Examples:
        >>> GriddedChord(Chord((0, 1, 1, 0)), ((0,0), (0,1), (2,1), (2,0))).get_bounding_indices(0, 1)
        0, 2, 2, 3
        >>> GriddedChord(Chord((0, 1, 1, 0)), ((0,0), (1,1), (1,1), (1,0))).get_bounding_indices(0, 1)
        0, 1, 1, 4
        >>> GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (1,0), (2,1))).get_bounding_indices(0, 2)
        0, 1, 3, 4
        """
        assert(col_source <= col_sink)

        points_left_col_source = list(self.get_points_left_col(col_source))
        points_right_col_source = list(self.get_points_right_col(col_source))
        # finds minimun index a source of a new chord in col_source can have, by checking what the largest index left of col_source
        min_index_source = max([idx for idx, _ in points_left_col_source] + [-1]) + 1
        # finds maximum index a source of a new chord in col_source, checking what the smallest index right col_source is.
        max_index_source = min([idx for idx, _ in points_right_col_source] + [len(self) * 2])

        points_left_col_sink = list(self.get_points_left_col(col_sink))
        points_right_col_sink = list(self.get_points_right_col(col_sink))

        # finds minimum index a sink of a new chord in row_sink can have, checking what the largest chord below row_sink is
        min_index_sink = max([idx for idx, _ in points_left_col_sink] + [min_index_source - 1]) + 1
        # finds maximum index a new sink chord can have in row_sink, checking what the smallest chord above row_sink is
        max_index_sink = min([idx for idx, _ in points_right_col_sink] + [len(self.patt)])

        return (min_index_source, max_index_source, min_index_sink, max_index_sink)
    
    # this was copied in for use in chord_placement
    def get_bounding_box(self, cell: Cell) -> Tuple[int, int, int, int]:
        """Determines the range of indices and values of the points in 
        the gridded chord diagram that can be found in the Cell cell."""
        row = list(self.get_points_row(cell[1]))
        col = list(self.get_points_col(cell[0]))
        if not row: # if there are no points in the selected row
            above = list(self.get_points_above_row(cell[1]))
            if not above:
                maxval = len(self)
            else:
                maxval = min(point[1] for point in above)
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
                maxdex = len(self.patt) 
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
        return (mindex, maxdex, minval, maxval)


    ### Modify instance of GriddedChord ###
    def insert_chord(self, row: int, col_source: int, col_sink: int) -> Iterator["GriddedChord"]:
        # sCN: insert_point changed to insert_chord, functionality changed appropriately
        """Insert a new point into cell of the gridded chord, such that the
        point is added to the underlying pattern with the position at the cell.
        Yields all gridded chords where the point has been mixed into the points
        in the cell.
        
        Examples:
        >>> GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 1), (0, 1))).insert_chord(0, 0, 1)
        [GriddedChord(Chord((0, 1, 2, 2, 0, 1)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 0, 2, 1)), ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 2, 1, 0)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 1, 2, 0)), ((0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 2, 1, 0)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0))),
        GriddedChord(Chord((0, 1, 2, 2, 0, 1)), ((0, 0), (0, 0), (0, 1), (1, 1), (1, 0), (1, 0)))]
        >>> GriddedChord(Chord((0, 1, 1, 0)), ((0, 0), (1, 0), (1, 1), (0, 1))).insert_chord(0, 1, 1)
        []
        """
        possible_chords = []
        min_index_source, max_index_source, min_index_sink, max_index_sink = self.get_bounding_indices(col_source, col_sink)
        
        # sAsk: There is probably a better way to do this
        for source in range(min_index_source, max_index_source + 1):
            for sink in range(max(source, min_index_sink), max_index_sink + 1):
                chord = self.insert_specific_chord(row, col_source, col_sink, source, sink)
                if chord != None:
                    possible_chords.append(chord)
        return possible_chords

    def insert_specific_chord(self, row: int, col_source: int, col_sink: int, source: int, sink: int) -> "GriddedChord":
        # sCN: insert_specific_point -> insert_specific_chord, changed numebr of parameters, adapted to chords appropriately
        #sToDo: code check for vaild gridded chord required
        """Insert a new chord in row, with source and sink in col_source and col_sink respecitively, at indices source and sink.
        """
        patt = self._chord
        patt = patt.insert(source, sink)
        positions = list(self._pos)
        positions.insert(source, (col_source, row))
        positions.insert(sink + 1, (col_sink, row)) # +1 needed to account for source being inserted first.
        try:
            gc =  GriddedChord(patt, positions)
            return gc
        except ValueError:
            return None
    
    def remove_chord(self, chord: int) -> "GriddedChord":
        #sCN: changed to work with chords
        """Remove the chord of the int given. Doesn't do anything if an incorrect chord is given.
        
        Examples:
        >>> GriddedChord(Chord((0,1,2,0,2,1)), ((0,0), (1,0), (1,1), (0,1), (1,2), (1,3))).remove_chord(0)
        GriddedChord(Chord((0,1,1,0), ((1,0), (1,1), (1,2), (1,3)))
        >>> GriddedChord(Chord((0,1,2,0,2,1)), ((0,0), (1,0), (1,1), (0,1), (1,2), (1,3))).remove_chord(3)
        GriddedChord(Chord((0,1,2,0,2,1)), ((0,0), (1,0), (1,1), (0,1), (1,2), (1,3)))
        >>> GriddedChord(Chord((0,1,2,0,2,1)), ((0,0), (1,0), (1,1), (0,1), (1,2), (1,3))).remove_chord(2)
        GriddedChord(Chord((0,1,0,1)), ((0,0), (1,0), (0,1), (1,3)))
        """
        patt = Chord.to_standard([item for item, pos in self if item != chord])
        pos = [pos for item, pos in self if item != chord]
        return type(self)(patt, pos)

    def remove_chord_idx(self, index) -> "GriddedChord":
    # sCN: used cleaner list notation
        """Returns the Gridded Chord with the chord at index removed.
        """
        avoid_chord = self._patt[index]
        new_chords = [chord for chord, _ in self if chord != avoid_chord]
        new_positions = [pos for chord, pos in self if chord != avoid_chord]

        new_patt = Chord.to_standard(new_chords)

        return GriddedChord(new_patt, new_positions)
    
    def remove_cells(self, cells: Iterable[Cell]) -> "GriddedChord":
    # sCN: old method was perm specific
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
    
    def apply_map(self, cell_mapping: Callable[[Cell], Cell]) -> "GriddedChord":
        """Map the coordinates to a new list of coordinates according to the
        cell_mapping given.
        
        Example:
        >>> def add_one_x(cell: Cell):
             return (cell[0] + 1, cell[1])

        >>> GriddedChord(Chord((0, 1, 0, 1)), ((0,0), (1,1), (0,1), (1,2))).apply_map(add_one_x)
        GriddedChord(Chord((0, 1, 0, 1)), ((1,0), (2,1), (1,1), (2,2)))"""
        return type(self)(self._chord, [cell_mapping(cell) for cell in self._pos])
    
    def apply_chord_map_to_cell(
        self, chord_mapping: Callable[[Chord], Chord], cell: Cell
    ) -> "GriddedChord":
    # skipped testing
    # looks about right, have not tested.
        """Apply a chord map to the subchord within a cell."""
        subchord = [val for val, pos in self if pos == cell]
        st = Chord.to_standard(subchord)
        back_map = dict(zip(st, subchord))
        new_subchord = (back_map[val] for val in chord_mapping(st))
        return type(self)(
            (next(new_subchord) if pos == cell else val for val, pos in self), self.pos
        )


    ### Miscellaneous/Not sure what is happening ###
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
    
    # sToDo: not sure what this is doing, need to look more at union find algorithm
    # sToDo: tests
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
        return [self.get_subchord_in_cells(comp) for comp in factor_cells]

    # this method does not seem relavent, also don't know what it is doing sAsk
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


    ### Data processing and printing ###
    @property
    def patt(self) -> Chord:
        return self._patt

    @property
    def pos(self) -> Tuple[Cell, ...]:
        return self._pos

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
            Chord(tuple(next(iterator) for _ in range(length))), ((next(iterator), next(iterator)) for _ in range(length))
        )

    def to_jsonable(self) -> dict:
        """Returns a dictionary object which is JSON serializable representing
        a GriddedChord."""
        return {"patt": self._patt, "pos": self._pos}

    @classmethod
    def from_json(cls, jsonstr: str) -> "GriddedChord":
    # not sure how to test (what is a JSON string)
        """Returns a GriddedChord object from JSON string."""
        jsondict = json.loads(jsonstr)
        return cls.from_dict(jsondict)

    @classmethod
    def from_dict(cls, jsondict: dict) -> "GriddedChord":
        """Returns a GriddedChord object from a dictionary loaded from a JSON
        serialized GriddedChord object."""
        return cls(Chord(jsondict["patt"]), map(tuple, jsondict["pos"]))  # type: ignore
    
    # tests to here
    def ascii_plot(self) -> str:
    # needs fixing for flipped chords
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
    

    ### Dunder methods ###
    def __len__(self) -> int:
        return len(self._chord)

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
