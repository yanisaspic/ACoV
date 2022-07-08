"""
API to facilitate the comparison between different embryos.
_______________

    + replace segmentation-specific time points with comparable minutes-post-fertilization
    + apply a voxelsize correction to account for experimental differences in the image acquisition
_______________

@ ASLOUDJ Yanis
07/08/2022
"""


import pandas as pd
import numpy as np
from sklearn.linear_model import RANSACRegressor
from os import listdir


def apply_to_dict_values(dic, fun, **kwargs):
    """
    # Description
    ---
    Applies a function to all the values (nested values included) of a dict.

    # Argument(s)
    ---
        `dic` (dict): a dictionary.
        `fun` (fun): a function to apply to all the dict entries.

    # Usage
    >>> dic, fun = {'a': 10, 'b': {'c': 30}}, lambda x: x+3
    >>> print(apply_to_dict_values(dic, fun))
    ... {'a': 13, 'b': {'c': 33}}
    """
    output = dic.copy()
    for k in output:
        if isinstance(output[k], dict):
            output[k] = apply_to_dict_values(output[k], fun, **kwargs)
        else:
            output[k] = fun(output[k], **kwargs)
    return output

def append_to_identical_dict(dic1, dic2, k='none'):
    """
    # Description
    ---
    Appends the values of a first dict into a second one with an identical key structure.

    # Argument(s)
    ---
        `dic1` (dict): the dictionary on which the data is appended.
        `dic2` (dict): the dictionary with the data to append.
        `k` (str): a variable called recursively.

    # Usage
    ---
    >>> dic1, dic2 = {'a': [], 'b': {'c': [3]}}, {'a': 10, 'b': {'c': 20}}
    >>> print(append_to_identical_dict(dic1, dic2))
    ... {'a': [10], 'b': {'c': [3, 20]}}
    """
    output = dic1.copy()
    for k in output:
        if isinstance(output[k], dict):
            append_to_identical_dict(output[k], dic2[k], k)
        else:
            output[k].append(dic2[k])
    return output

def concatenate_embryo_resolution_data(xlsx_folder):
    """
    # Description
    ---
    Concatenates the embryo-resolution data of all the .xlsx individual files found in the corresponding folder.
    A new embryo column is added to the data frame: it corresponds to the name of the embryo data source.

    # Usage
    ---
    >>> all_embryo_resolution_data = concatenate_embryo_resolution_data('xlsx')
    >>> print(all_embryo_resolution_data)
    ...    tp  object    volume  total_surface  cell_count     embryo
        0    0  embryo  67400259   9.063063e+05          46  Astec-pm3
        1    1  embryo  67789699   8.963664e+05          47  Astec-pm3
        2    2  embryo  67878724   9.167023e+05          53  Astec-pm3
        3    3  embryo  68163715   9.325861e+05          60  Astec-pm3
        4    4  embryo  68172491   9.078513e+05          62  Astec-pm3
        ..  ..     ...       ...            ...         ...        ...
        85  85  embryo  60935962   1.040680e+06         299  Astec-pm5
        86  86  embryo  60729570   1.037320e+06         305  Astec-pm5
        87  87  embryo  60447756   1.041947e+06         311  Astec-pm5
        88  88  embryo  60374948   1.038938e+06         322  Astec-pm5
        89  89  embryo  60304297   1.032776e+06         329  Astec-pm5
    """
    all_embryo_resolution_data = []
    for filename in listdir(xlsx_folder):
        embryo_name = filename.split('.')[0]
        individual_embryo_resolution_data = pd.read_excel(f'{xlsx_folder}/{filename}', engine='openpyxl', sheet_name='embryo')
        individual_embryo_resolution_data['embryo'] = embryo_name
        all_embryo_resolution_data.append(individual_embryo_resolution_data)
    return pd.concat(all_embryo_resolution_data).reset_index(drop=True)

def get_volume_coefficients(all_embryo_resolution_data):
    """
    # Description
    ---
    

    # Argument(s)
    ---
        `all_embryo_resolution_data` (pandas df): data structure generated using concatenate_embryo_resolution_data().
    """

print(concatenate_embryo_resolution_data('xlsx'))