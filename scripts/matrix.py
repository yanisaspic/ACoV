"""
API to convert an .xlsx file into a matrix of cell composition during the embryogenesis.
_______________

    + binary cell presence-absence matrices
    + tissue size expressed in cells matrices
_______________

@ ASLOUDJ Yanis
07/11/2022
"""


from scripts.alignment import concatenate_unique_resolution_data

import pandas as pd


def fill_dataframe_columns(df):
    """
    # Description
    ---
    Adds intermediate columns to an existing dataframe with integer numbers as column names.
    The values of the new columns are identical to the column before them.

    # Argument(s)
    ---
        `df` (pandas df): a dataframe with integer numbers as columns names.

    # Usage
    ---
    >>> data = [{1: 10, 2: 12, 4: 28, 7:31}, {1: 14, 2: 10, 4: 6, 7:1}]
    >>> df = pd.DataFrame.from_records(data)
    >>> print(df)
    ...    1   2   4   7
        0  10  12  28  31
        1  14  10   6   1
    >>> print(fill_dataframe_columns(df))
    ...    1   2   3   4   5   6   7
        0  10  12  12  28  28  28  31
        1  14  10  10   6   6   6   1
    """
    tmp = df.copy()
    left_col = tmp.columns[0]

    for col in tmp:
        right_col = col
        col_step = right_col - left_col

        while col_step>1:
            tmp[left_col+1] = tmp[left_col]
            left_col += 1
            col_step = right_col - left_col
        left_col = right_col

    return tmp.sort_index(axis=1)

def fill_dataframe_rows(df, row_names):
    """
    # Description
    ---
    Adds new named rows to an existing dataframe.
    The rows values are 0. They are sorted ascendingly.

    # Argument(s)
    ---
        `df` (pandas df): a dataframe.
        `row_names` (list): names of the new rows to add.

    # Usage
    ---
    >>> data = [{1: 10, 2: 12, 4: 28, 7:31}, {1: 14, 2: 10, 4: 6, 7:1}]
    >>> df = pd.DataFrame.from_records(data)
    >>> print(df)
    ...    1   2   4   7
        0  10  12  28  31
        1  14  10   6   1
    >>> print(fill_rows(df, [12, 8, 6]))
    ...     1   2   4   7
        0   10  12  28  31
        1   14  10   6   1
        6    0   0   0   0
        8    0   0   0   0
        12   0   0   0   0
    """
    tmp = df.copy()
    tmp = tmp.T
    for index in row_names:
        tmp[index] = 0
    tmp = tmp.T
    return tmp.sort_index(axis=0)

def get_cell_count_matrices(sub_resolution_data):
    """
    # Description
    ---
    Returns a dict containing matrices with objects as rows, a time axis proxy as column, and the cell_counts as values.
    The dict is split into two time axes subdicts: 'embryo_cell_count' and 'minutes_post_fertilization'.
    Each subdict associates an embryo name to a matrix.

    # Argument(s)
    ---
        `sub_resolution_data` (pandas df): 'tissue' or 'cell' sheet in a .xlsx file.

    # Usage
    ---
    >>> cell_resolution_data = concatenate_unique_resolution_data('xlsx', 'cell')
    >>> matrices = get_cell_count_matrices(cell_resolution_data)
    >>> print(matrices['embryo_cell_count']['Astec-pm8'])
    ... embryo_cell_count  64   65   66   ...  395  396  397                                                       
        A10.25*            0.0  0.0  0.0 ...  1.0  1.0  1.0
        A10.25_            0.0  0.0  0.0 ...  1.0  1.0  1.0
        A10.26*            0.0  0.0  0.0 ...  1.0  1.0  1.0
        A10.26_            0.0  0.0  0.0 ...  1.0  1.0  1.0
        A10.27*            0.0  0.0  0.0 ...  1.0  1.0  1.0
    """
    matrices = {'embryo_cell_count': {}, 'minutes_post_fertilization': {}}
    all_objects = set(sub_resolution_data.object.unique())

    for embryo_name in sub_resolution_data.embryo.unique():
        individual_sub_resolution_data = sub_resolution_data[sub_resolution_data.embryo==embryo_name]

        for time_axis in matrices.keys():
            matrix = individual_sub_resolution_data.pivot_table(index='object', columns=time_axis, values='cell_count')
            missing_objects = all_objects.difference(matrix.index)

            matrix = fill_dataframe_rows(matrix, missing_objects)
            matrix = fill_dataframe_columns(matrix)
            matrices[time_axis][embryo_name] = matrix.fillna(0)
            matrix.index.name = None

    return matrices

def save_cell_count_matrices(matrices_dict, filename):
    """
    # Description
    ---
    Saves a subdict of cell_count matrices into an .xlsx file where each sheet corresponds to the matrix of an embryo.

    # Argument(s)
    ---
        `matrices_dict` (dict): embryo names as keys and cell count matrices as values.
        `filename` (str): name of the resulting .xlsx file.
    """
    with pd.ExcelWriter(f'{filename}.xlsx') as writer:
        for embryo_name, matrix in matrices_dict.items():
            matrix.to_excel(writer, sheet_name=embryo_name)