import string
from typing import Tuple, List, Optional, Union, Generator
import csv
from pathlib import Path
from deprecated import deprecated


class Plate:
    """
    A plate consists of a rectangular grid of wells and an index for rows and columns.
    The row indices are uppercase characters starting from A.
    The column indices are integers starting at 1.
    Each well is defined as a string row + column, e.g. 'A1'
    --------------------
    Parameters:
        n_rows: Number of rows
        n_cols: Number of columns
        max_vol: Maximum volume per well in nL
        dead_vol: Dead volume per well in nL
    --------------------
    The state of a plate instance is fully defined by the compounds and volume every well holds.
    """
    def __init__(self, n_rows: int, n_cols: int, max_vol: int, dead_vol: int):
        if dead_vol > max_vol:
            raise ValueError('Dead volume cannot exceed maximum volume')
        if n_rows < 1 or n_cols < 1:
            raise ValueError('Invalid plate format')
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.max_vol = max_vol  # max max_vol per well in nL
        self.dead_vol = dead_vol  # dead volume per well in nL
        self._compounds: List[List[List[str]]] = [[[] for _ in range(n_cols)] for _ in range(n_rows)]
        self._volume: List[List[int]] = [[0 for _ in range(n_cols)] for _ in range(n_rows)]

    def __str__(self):
        """Print an overview of the Plate"""
        printout = '  '
        printout += ''.join([f'{i + 1}'.ljust(32) for i in range(self.n_cols)])
        printout += '\n'
        for row_n, (row_c, row_v) in enumerate(zip(self._compounds, self._volume)):
            printout += chr(65+row_n)
            for elem_c, elem_v in zip(row_c, row_v):
                printout += f' ({elem_c}, {elem_v})'.ljust(32)
            printout += '\n'
        return printout

    def __sanitize_inputs(self, row: int, col: int, compounds: Union[str, List[str]], vol: int):
        """
        Check if user-entered values make sense for the limitations of the current plate.
        Dead volume limitations are not checked here, as some methods (e.g. Plate.empty_well() ar supposed to violate
        this criterion. Check the dead volume limitation in public methods where necessary before passing the
        user-entered values to a private method.
        """
        if row > self.n_rows or row < 0:
            raise IndexError(f'Row index is out-of-bounds for plate with {self.n_rows} rows')
        if col > self.n_cols or col < 0:
            raise IndexError(f'Column index is out-of-bounds for plate with {self.n_cols} columns')
        if vol < 0:
            raise ValueError('Volume cannot be negative')
        elif vol > self.max_vol:
            raise ValueError('Volume cannot exceed maximum volume')
        if vol > 0 and len(compounds) < 1:
            raise ValueError('Cannot assign non-zero volume without assigning compound')
        if type(compounds) is list:
            compounds = [c.strip() for c in compounds]
        else:
            compounds = compounds.strip()
        if compounds == '' or compounds == ['']:
            """this fixes a bug where loading/saving would of a plate would change it and provides a standard for how 
            an empty well should look like"""
            compounds = []
        return row, col, compounds, vol

    def __indices(self) -> Tuple[List, List]:
        """Return the row and column indices. For use by other private methods."""
        return list(range(self.n_rows)), list(range(self.n_cols))

    def __get_well(self, row: int, col: int):
        """Return list of compounds and volumes in one well. Base method to inspect the plate"""
        return self._compounds[row][col], self._volume[row][col]

    def __set_well(self, row: int, col: int, compounds: List[str], vol: int):
        """Set the compounds and volumes in one well. Base method to modify the plate"""
        row, col, compounds, vol = self.__sanitize_inputs(row, col, compounds, vol)
        self._compounds[row][col] = compounds
        self._volume[row][col] = vol

    def __get_span(self, rows: list, cols: list) -> Tuple[List, List]:
        """
        For a list of rows and columns, return the compounds and volumes in all wells spanned by that rows and columns.
        The returned tuple contains two lists of shape (n_rows x n_columns)
        """
        return [[self._compounds[row][col] for col in cols] for row in rows], [
            [self._volume[row][col] for col in cols] for row in rows]

    def __set_span(self, rows: list, cols: list, compounds: list, volumes: list):
        """
        For a list of rows and columns, set the compounds and volumes in all wells spanned by that rows and columns.
        The passed lists must be of shape (n_rows x n_columns)
        """
        # validate input
        try:
            assert len(compounds) == len(rows)
            assert len(volumes) == len(rows)
            for row in compounds:
                assert len(row) == len(cols)
            for row in volumes:
                assert len(row) == len(cols)
        except AssertionError:
            raise ValueError('Input shape not appropriate')
        # set all wells in the span
        for r, row in enumerate(rows):
            for c, col in enumerate(cols):
                self.__set_well(row, col, compounds[r][c], volumes[r][
                    c])  # We need the enumeration because row and col refer to the entire plate, but compounds and volumes only refer to the span, so using row and col raises IndexErrors

    # TODO deprecate methods that are redundant with __get_span and __set_span

    @deprecated(reason='Use __get_span() instead')
    def __get_row(self, row: int) -> Tuple[List, List]:
        row = [row, ]
        cmp, vol = self.__get_span(row, self.__indices()[1])
        return cmp[0], vol[0]

    @deprecated(reason='Use __get_span() instead')
    def __get_column(self, col: int):
        return [self._compounds[row][col] for row in range(self.n_rows)], [
            self._volume[row][col] for row in range(self.n_rows)]

    def __is_filled(self, row: int, col: int) -> bool:
        vol = self._volume[row][col]
        if vol == 0:
            return False
        else:
            return True

    @deprecated(reason='Use __set_span() instead')
    def __set_row(self, row: int, compounds: List[list], vol: List[int]):
        for col in range(self.n_cols):
            self.__set_well(row, col, compounds[col], vol[col])

    @deprecated(reason='Use __set_span() instead')
    def __set_column(self, col: int, compounds: List[list], vol: List[int]):
        for row in range(self.n_rows):
            self.__set_well(row, col, compounds[row], vol[row])

    @deprecated(reason='Use __set_span() instead')
    def __set_plate(self, compounds: list, vol: int):
        for row in range(self.n_rows):
            for col in range(self.n_cols):
                self.__set_well(row, col, compounds, vol)

    @staticmethod
    def __to_index(pos_str: str) -> Tuple[int, int]:
        """Convert a location string of the form [A-Z][0-9]{n} to indices for row and column """
        row_str, col_str = pos_str[0], pos_str[1:]
        row_idx = string.ascii_uppercase.index(row_str)
        col_idx = int(col_str) - 1
        return row_idx, col_idx

    @staticmethod
    def __from_index(row_idx: int, col_idx: int) -> str:
        """Reverse method of self.__to_index()"""
        row_str = string.ascii_uppercase[row_idx]
        col_str = str(col_idx + 1)
        return ''.join([row_str, col_str])

    def shape(self) -> Tuple[int, int]:
        """Return dimension of the plate"""
        return self.n_rows, self.n_cols

    def rows(self) -> Tuple[str, ...]:
        """Return the row identifiers of the plate"""
        return tuple(self.__from_index(row_idx, 0)[0] for row_idx in range(self.n_rows))

    def columns(self) -> Tuple[str, ...]:
        """Return the column identifiers of the plate"""
        return tuple(self.__from_index(0, col_idx)[1:] for col_idx in range(self.n_cols))

    def wells(self) -> Tuple[str, ...]:
        """Return the well identifiers of the plate"""
        return tuple(self.__from_index(row_idx, col_idx)
                     for row_idx in range(self.n_rows) for col_idx in range(self.n_cols))

    def fill_well(self, pos: str, compound: Union[str, list], vol: int):
        """Add single compound to well"""
        row, col = self.__to_index(pos)
        cmp_cur, vol_cur = self.__get_well(row, col)
        if type(compound) is str:
            compound = [compound]  # turn into list
        self.__set_well(row, col, cmp_cur + compound, vol_cur + vol)

    def set_well(self, pos: str, compounds: List[str], vol: int):
        """Set well content, overwriting any previous content"""
        row, col = self.__to_index(pos)
        self.__set_well(row, col, compounds, vol)

    def consume_well(self, pos: str, vol: int, force=False):
        """
        Remove a given volume from well
        :param pos: Well to remove aliquot from
        :param vol: Volume to be removed
        :param force: default False. if True, the dead volume limitation of the well can be violated
        """
        row, col = self.__to_index(pos)
        cmp_cur, vol_cur = self.__get_well(row, col)
        if vol_cur < self.dead_vol and force is False:
            raise ValueError('Cannot consume well below dead volume level')
        self.__set_well(row, col, cmp_cur, vol_cur - vol)

    def empty_well(self, pos):
        """Set well volume to zero and delete compounds"""
        row, col = self.__to_index(pos)
        self.__set_well(row, col, [], 0)
    # TODO fill_span(start_well, end_well, cmp, vol) could superseed fill_row, fill_column and fill_block

    def fill_span(self, well_start: str, well_end: str, compound: str, vol: int):
        """Add single compound to all wells in the rectangle spanned by well_start and well_end"""
        row_start, col_start = self.__to_index(well_start)
        row_end, col_end = self.__to_index(well_end)
        rows = list(range(row_start, row_end + 1))
        cols = list(range(col_start, col_end + 1))
        cur_cmp, cur_vol = self.__get_span(rows, cols)
        new_cmp = [[well + [compound] for well in row] for row in cur_cmp]
        new_vol = [[well + vol for well in row] for row in cur_vol]
        self.__set_span(rows, cols, new_cmp, new_vol)

    def empty_span(self, well_start: str, well_end: str):
        """Set volume to 0 and delete compounds in all wells in the rectangle spanned by well_start and well_end"""
        row_start, col_start = self.__to_index(well_start)
        row_end, col_end = self.__to_index(well_end)
        rows = list(range(row_start, row_end + 1))
        cols = list(range(col_start, col_end + 1))
        cur_cmp, cur_vol = self.__get_span(rows, cols)
        new_cmp = [[[] for well in row] for row in cur_cmp]
        new_vol = [[0 for well in row] for row in cur_vol]
        self.__set_span(rows, cols, new_cmp, new_vol)

    def fill_column(self, col: str, compound: str, vol: int):  # TODO rewrite to use span logic
        """Add single compound to all wells in a column"""
        col = self.rows()[0] + col  # auxiliary row, no influence
        col = self.__to_index(col)[1]
        cmp_cur_list, vol_cur_list = self.__get_column(col)
        self.__set_column(col, [cmp + [compound] for cmp in cmp_cur_list], [v + vol for v in vol_cur_list])

    def fill_block(self, rows: Tuple[str], cols: Tuple[str], compound: str, vol: int):  # TODO rewrite to use span logic
        """Add single compound to all wells in a block that spans one or more rows and one or more columns"""
        for row_str in rows:
            for col_str in cols:
                self.fill_well(row_str + col_str, compound, vol)

    def fill_plate(self, compound: str, vol: int):
        """Fill all wells of a plate with identical compound and volume"""
        self.__set_plate([compound], vol)

    def consume_plate(self, vol: int, force=False):
        """Remove a given volume from every well of the plate"""
        for row in range(self.n_rows):
            for col in range(self.n_cols):
                cmp_cur, vol_cur = self.__get_well(row, col)
                if vol_cur < self.dead_vol and force is False:
                    raise ValueError('Cannot consume well below dead volume level')
                self.__set_well(row, col, cmp_cur, vol_cur - vol)

    def empty_plate(self):
        """Set all well volumes to zero and delete compounds."""
        self.__set_plate([], 0)

    # def write_plate(self, compounds, vols):
    #     """Set compounds and volumes on entire plate. This overwrites previous contents."""
    #     self.__set_span(self.__indices()[0],
    #                     self.__indices()[1],
    #                     )

    def well(self, pos) -> tuple:
        """Return the compounds and volumes in a well"""
        row, col = self.__to_index(pos)
        return self.__get_well(row, col)

    def compounds(self, pos) -> list:
        """Return the list of compounds in a well"""
        row, col = self.__to_index(pos)
        return self.__get_well(row, col)[0]

    def volume(self, pos) -> int:
        """Return the total volume in a well"""
        row, col = self.__to_index(pos)
        return self.__get_well(row, col)[1]

    def iterate_wells(self) -> Generator[tuple, None, None]:
        """Yield a tuple (compounds, volume) for every well in the plate"""
        for well in self.wells():
            yield self.well(well)

    def row(self, row: str) -> Tuple[List[List[str]], List[int]]:
        """Return compounds and volumes in a row"""
        compounds, vols = [], []
        for col in self.columns():
            cmp, vol = self.well(row + col)
            compounds.append(cmp)
            vols.append(vol)
        return compounds, vols

    def column(self, col: str) -> Tuple[List[List[str]], List[int]]:
        """Return compounds and volumes in a column"""
        compounds, vols = [], []
        for row in self.rows():
            cmp, vol = self.well(row + col)
            compounds.append(cmp)
            vols.append(vol)
        return compounds, vols

    def free(self) -> Optional[str]:
        """Return the position of the next free well. Order is left to right, top to bottom."""
        for well in self.wells():
            if self.volume(well) == 0:
                return well
        return None

    def to_csv(self, file, save_volumes=False):
        """Save the components in the plate into a csv file. Volumes are not saved by default."""
        if isinstance(file, Path): # convert any Path to str
            file = str(file.resolve())

        with open(file, 'w') as csv_file:
            writer = csv.writer(csv_file)
            # write header
            header = ['', ]
            header += self.columns()
            writer.writerow(header)
            for row_str in self.rows():
                text = [f'{row_str}', ] + [', '.join(elem) for elem in self.row(row_str)[0]]
                writer.writerow(text)
        if save_volumes:
            vol_file = file.strip('.csv') + '_volumes.csv'
            with open(vol_file, 'w') as csv_file:
                writer = csv.writer(csv_file)
                # write header
                header = ['', ]
                header += self.columns()
                writer.writerow(header)
                for row_str in self.rows():
                    text = [f'{row_str}', ] + [str(elem) for elem in self.row(row_str)[1]]
                    writer.writerow(text)

    def from_csv(self, file, vol=None):
        """
        Load a plate from csv file. The csv file must be structured analogous to the output of
        plate.to_csv(save_volumes=True).
        Currently, 96 and 384 well plates are supported.
        :param file: str or Path-like Path to the csv file. A corresponding volume file must be in the same directory
        The shape must correspond to the shape of self. Alternatively, if the volume is identical for all wells,
        a optional param vol may be given instead.
        :param vol: int, optional parameter. If vol is given, volume files will be ignored and all wells of the
        target plate will be initialized to this volume
        WARNING: Compounds / volumes in self will be overwritten.
        """
        # TODO make this infer whether we need 96 or 384 (or1536?) from supplied data
        # TODO it would be better not to have this as a class method, so that it can return a Plate. (like pd.read_csv())
        if isinstance(file, Path):  # convert any Path to str
            file = str(file.resolve())

        vol_file = file.strip('.csv') + '_volumes.csv'
        """parse files"""
        parsed = []
        with open(file, 'r') as csv_file:  # compound file
            reader = csv.reader(csv_file)
            content = [row for row in reader]
        col_str = content.pop(0)  # the first line must be column names
        col_str.pop(0)  # the first item (emtpy string) is irrelevant
        row_str = [row.pop(0) for row in content]
        parsed.append([row_str, col_str, content])

        if vol is None:  # standard case: a volume file was supplied but no optional parameter vol
            with open(vol_file, 'r') as csv_file:  # same as above for vol file
                reader = csv.reader(csv_file)
                content = [row for row in reader]
            col_str = content.pop(0)  # the first line must be column names
            col_str.pop(0)  # the first item (emtpy string) is irrelevant
            row_str = [row.pop(0) for row in content]
            parsed.append([row_str, col_str, content])
        else:  # optional parameter vol was supplied: volume file will be ignored and does not need to exist
            # Map the volume to all filled wells (use 0 for empty wells)
            static_volumes = [[vol if elem else 0 for elem in row] for row in content]
            parsed.append([row_str, col_str, static_volumes])

        """control if inputs are valid"""
        if parsed[0][0] != parsed[1][0] or parsed[0][1] != parsed[1][
            1]:  # equal dimensions for compound and volume files
            raise ValueError(f'Plate shape in {file} does not fit plate shape in {vol_file}')
        else:
            row_str, col_str, compounds, volumes = parsed[0][0], parsed[0][1], parsed[0][2], parsed[1][2]

        # rows are alphabetic values from 0 (A) to 8*n-1
        valid_rows = [[string.ascii_uppercase[i] for i in range(0, 8)],
                      [string.ascii_uppercase[i] for i in range(0, 16)],
                      ]
        # cols are int from 1 to 12*n
        valid_cols = [[str(i + 1) for i in range(0, 12)],
                      [str(i + 1) for i in range(0, 24)],
                      ]

        if tuple(row_str) != self.rows():
            raise ValueError('Invalid row')
        if tuple(col_str) != self.columns():
            raise ValueError('Invalid column')
        compounds = [[s.split(sep=',') for s in row] for row in compounds]
        volumes = [[int(i) for i in row] for row in volumes]  # cast to int
        """now that the plate format is validated, set the new plate"""
        self.__set_span(self.__indices()[0],
                        self.__indices()[1],
                        compounds,
                        volumes
                        )

    def to_dict(self, include_volumes=False):
        """
        Convert plate to dictionary. Well identifiers are used as keys, a list of contents is used as values.
        :param include_volumes: bool, toggles whether volumes should be included in the dictionary
        :return: dict

        Example returns:
            include_volumes = False:
                > {'A1': ['water', 'EtOH'], 'A2': []}
            include_volumes = True
                > {'A1': (['water', 'EtOH'], 1000), 'A2': ([], 0)}
        """
        dic = {}
        if include_volumes is True:
            for well in self.wells():
                dic[well] = self.well(well)
        else:
            for well in self.wells():
                dic[well] = self.compounds(well)
        return dic

class Plate96(Plate):
    def __init__(self, max_vol: int, dead_vol: int):
        super().__init__(8, 12, max_vol, dead_vol)


class Plate384(Plate):
    def __init__(self, max_vol: int, dead_vol: int):
        super().__init__(16, 24, max_vol, dead_vol)


class Plate384Echo(Plate384):
    def __init__(self):
        super().__init__(max_vol=12000, dead_vol=2500)


if __name__ == '__main__':
    # debugging commands go here
    pass







