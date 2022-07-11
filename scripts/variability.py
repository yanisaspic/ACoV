"""
API to quantify the inter-individual variability during the embryogenesis.
The default parameters are described in the script utils.variability.py.
_______________

    + min-max normalized evenness indices
    + mean pairwise distances normalized by the embryo cell count
    + miscellaneous methods to verify the robustness
_______________

@ ASLOUDJ Yanis
07/11/2022
"""


import pandas as pd
import numpy as np
from mlxtend.preprocessing import TransactionEncoder


def one_hot_encode(dataset):
    """
    # Description
    ---
    Applies a one-hot encoding on a dataset of lists and returns the result as a dataframe.

    # Argument(s)
    ---
        `dataset` (list of lists): each sublist has a label.

    # Usage
    ---
    >>> dataset = [ ['B', 'K', 'L'], ['L', 'K'], ['B', 'L'] ]
    >>> print(one_hot_encode(dataset))
    ...   B  K  L
        0  1  1  1
        1  0  1  1
        2  1  0  1
    """
    te = TransactionEncoder()
    te_ary = te.fit(dataset).transform(dataset)
    return pd.DataFrame(te_ary.astype('int'), columns=te.columns_)

def frame_embryo_coexistence(matrices_dict):
    """
    # Definition
    ---
    Returns a binary dataframe with embryos as rows and a time axis proxy as column.
    If measurements exist for an embryo at a given time, the corresponding cell has a value of 1.

    # Argument(s)
    ---
        `matrices_dict` (dict): embryo names as keys and cell count matrices as values.

    # Usage
    ---
    >>> data = pd.read_excel('cell_ecc.xlsx', engine='openpyxl', sheet_name=None, index_col=0)
    >>> print(frame_embryo_coexistence(data))
    ...        35   36   37   ...  695  696  697
    Astec-pm3    0    0    0  ...    0    0    0
    Astec-pm4    1    1    1  ...    0    0    0
    Astec-pm8    0    0    0  ...    0    0    0
    Astec-pm7    0    0    0  ...    0    0    0
    Louis        1    1    1  ...    1    1    1
    Astec-pm1    0    0    0  ...    1    1    1
    Astec-pm9    0    0    0  ...    0    0    0
    Astec-pm5    0    0    0  ...    0    0    0
    """
    coexistence = one_hot_encode(list(map( lambda matrix: list(matrix.columns), matrices_dict.values() )))
    coexistence.index = matrices_dict.keys()
    return coexistence.loc[:, coexistence.sum()>=2]

def get_distances_row(matrices_dict, pairwise_distance_method='cityblock', aggregate_distance_method='mean'):
    """
    # Description
    ---
    Returns a pandas Series corresponding to 
    """

def get_abundance_matrix(cell_matrices_dict):
    """
    # Description
    ---
    Returns an abundance matrix from a dict of binary cell-presence matrices.
    In the matrix, the rows are cells, the column is a time axis proxy and the value is the number of individuals with a given cell at a given time.
    The abundance is required to compute the evenness.

    # Argument(s)
    ---
        `cell_matrices_dict` (dict): embryo names as keys and cell count matrices as values at the cell-resolution.

    # Usage
    ---
    >>> data = pd.read_excel('cell_ecc.xlsx', engine='openpyxl', sheet_name=None, index_col=0)
    >>> print(get_abundance_matrix(data))
    ...
    """
    abundance = {}
    embryo_coexistence = frame_embryo_coexistence(cell_matrices_dict)
    time_points = embryo_coexistence.columns

    for binary_matrix in cell_matrices_dict.values():
        print(abundance.get(35, np.array([0] * len(time_points))))
        for time in time_points:
            abundance[time] = abundance.get(time, np.array([0] * len(time_points))) + np.array(binary_matrix[time])
    return pd.DataFrame.from_records(abundance)

data = pd.read_excel('cell_ecc.xlsx', engine='openpyxl', sheet_name=None, index_col=0)
print(get_abundance_matrix(data))

def get_evenness_row(abundance_matrix, ):
    """
    # Description
    ---
    
    """

def min_max_normalize_dataframe(df):
    """
    # Description
    ---
    Min-max normalize every row of a data frame.
    """

data = pd.read_excel('cell_ecc.xlsx', engine='openpyxl', sheet_name=None)
print(get_abundance_matrix(data))

data = pd.read_excel('cell_ecc.xlsx', engine='openpyxl', sheet_name=None)
print(frame_embryo_coexistence(data))