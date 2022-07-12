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


from utils.variability_settings import *

import pandas as pd
import scipy.spatial.distance as ssd
from itertools import combinations
from skbio.diversity import alpha_diversity
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

def get_distance_row(matrices_dict, pairwise_distance_method='cityblock', aggregate_method='mean'):
    """
    # Description
    ---
    Returns a dataframe corresponding to the aggregated pairwise distances of embryo cell composition.
    The dataframe has a single row.

    # Argument(s)
    ---
        `matrices_dict` (dict): embryo names as keys and cell count matrices as values.
        `pairwise_distance_method` (str): a valid ssd.pdist distance calculation method.
        `aggregate_method` (str): one of 'mean', 'min' or 'max'.

    # Usage
    ---
    >>> data = pd.read_excel('cell_ecc.xlsx', engine='openpyxl', sheet_name=None, index_col=0)
    >>> print(get_distance_row(data))
    ...
    """
    pairwise_distance_rows = []
    embryo_coexistence = frame_embryo_coexistence(matrices_dict)
    time_points = embryo_coexistence.columns

    for time in time_points:
        embryo_pool_at_time = embryo_coexistence[embryo_coexistence[time]==1].index

        for embryo_pair in combinations(embryo_pool_at_time, 2):
            u, v = matrices_dict[embryo_pair[0]][time], matrices_dict[embryo_pair[1]][time]
            row = {'time': time, 'distance': ssd.pdist((u, v), pairwise_distance_method), 'method': pairwise_distance_method}
            pairwise_distance_rows.append(row)

    pairwise_distance_df = pd.DataFrame.from_records(pairwise_distance_rows)
    distance_row = pairwise_distance_df.pivot_table(values='distance', index='method', columns='time', aggfunc=aggregate_method)
    return distance_row

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
    >>> print(get_abundance_matrix(data).head())
    ...
    """
    abundance = {}
    embryo_coexistence = frame_embryo_coexistence(cell_matrices_dict)
    time_points = embryo_coexistence.columns

    for binary_matrix in cell_matrices_dict.values():
        for time in time_points:
            try:
                abundance[time] += binary_matrix[time]
            except KeyError:
                abundance[time] = binary_matrix[time]

    return pd.DataFrame.from_records(abundance)

def get_evenness_row(cell_matrices_dict, evenness_index='simpson_e'):
    """
    # Description
    ---
    Returns a dataframe corresponding to the evenness of embryo cell composition.
    The evenness value ranges from 0 to 1.
    A high evenness value indicates that the embryos share a similar cell composition.
    The dataframe has a single row.

    # Argument(s)
    ---
        `cell_matrices_dict` (dict): embryo names as keys and cell count matrices as values at the cell-resolution.
        `evenness_index` (str): a valid skbio evenness index.

    # Usage
    ---
    >>> data = pd.read_excel('cell_ecc.xlsx', engine='openpyxl', sheet_name=None, index_col=0)
    >>> print(get_evenness_row(data))
    ...
    """
    abundance = get_abundance_matrix(cell_matrices_dict)
    evenness_rows = []
    for time in abundance.columns:
        evenness_rows.append({'time': time, 'evenness': alpha_diversity(evenness_index, abundance[time]), 'method': evenness_index})
    evenness_df = pd.DataFrame.from_records(evenness_rows)
    return evenness_df.pivot_table(values='evenness', index='method', columns='time')

def min_max_normalize(df):
    """
    # Description
    ---
    Min-max normalize every row of a data frame.
    X = (X - X_min) / (X_max - X_min)

    # Argument(s)
    ---
        `df` (pandas df): a dataframe with numerical values only.

    # Usage
    ---
    >>> dataset = [ [3, 5, 2], [4, 8, -5] ]
    >>> df = pd.DataFrame.from_records(dataset)
    >>> print(min_max_normalize(df))
    ...          0    1    2
        0  0.333333  1.0  0.0
        1  0.692308  1.0  0.0
    """
    df = df.T
    for col in df:
        df[col] = (df[col] - df[col].min()) / (df[col].max() - df[col].min())
    return df.T

def get_variability_dataframe(cell_matrices_dict, variability_fun, variability_methods_pool):
    """
    # Description
    ---
    Returns a min-max normalized matrix of cell composition variability.
    The variability can be computed using distances or evenness.
    A row is a quantification method and a column is a time in the embryo development.

    # Argument(s)
    ---
        `cell_matrices_dict` (dict): embryo names as keys and cell-resolution cell count matrices as values.
        `variability_fun` (fun): get_evenness_row or get_distance_row
        `variability_methods_pool` (list of str): valid quantification methods (i.e. ssd.pdist distances or skbio evenness indices).

    # Usage
    ---
    >>> data = pd.read_excel('cell_ecc.xlsx', engine='openpyxl', sheet_name=None, index_col=0)
    >>> print(get_variability_dataframe(data, get_distance_row, pairwise_distance_pool).head())
    ...
    """
    variability_rows = []
    for method in variability_methods_pool:
        variability_rows.append(variability_fun(cell_matrices_dict, method))
    variability_df = pd.concat(variability_rows)
    return min_max_normalize(variability_df)

def get_tissue_subgroup_resolution_distances(tissue_matrices_dict, tissue_subgroups=germ_lines):
    """
    # Description
    ---
    Returns a matrix quantifying the variability due to specific subtissues (e.g. germ lines).
    A row is a subtissue and a column is a time in the embryo development.

    # Argument(s)
    ---
        `tissue_matrices_dict` (dict): embryo names as keys and tissue-resolution cell count matrices as values.
        `tissue_subgroups` (dict): tissues subgroups names as keys and list of included tissues as values.

    # Usage
    ---
    >>> data = pd.read_excel('tissue_ecc.xlsx', engine='openpyxl', sheet_name=None, index_col=0)
    >>> print(get_tissue_subgroup_resolution_distances(data, germ_lines))
    ...
    """
    distances_rows = []
    for tissues in tissue_subgroups.values():
        tmp_dict = tissue_matrices_dict.copy()
        for embryo_name, matrix in tmp_dict.items():
            tmp_dict[embryo_name] = matrix.loc[tissues]
        distances_rows.append(get_distance_row(tmp_dict))

    tissue_subgroups_distances = pd.concat(distances_rows)
    tissue_subgroups_distances.index = tissue_subgroups.keys()
    return tissue_subgroups_distances



data = pd.read_excel('tissue_ecc.xlsx', engine='openpyxl', sheet_name=None, index_col=0)
print(get_tissue_subgroup_resolution_distances(data, germ_lines))
