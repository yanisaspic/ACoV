"""
API to facilitate the comparison between different embryos.
The references variables for the alignments are described in the script utils.alignment_settings.py.
_______________

    + replace segmentation-specific time points with comparable minutes-post-fertilization
    + apply a voxelsize correction to account for experimental differences in the image acquisition
_______________

@ ASLOUDJ Yanis
07/08/2022
"""


from scripts.utils.alignment_settings import *

import pandas as pd
import numpy as np
from sklearn.linear_model import RANSACRegressor
from os import listdir


def concatenate_embryo_resolution_data(xlsx_folder):
    """
    # Description
    ---
    Concatenates the embryo-resolution data of all the .xlsx individual files found in the corresponding folder.
    A new embryo column is added to the data frame: it corresponds to the name of the embryo data source.

    # Argument(s)
    ---
        `xlsx_folder` (str): name of the folder containing the parsed segmented embryos.

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
    embryo_resolution_data = []
    for filename in listdir(xlsx_folder):
        embryo_name = filename.split('.')[0]
        individual_embryo_resolution_data = pd.read_excel(f'{xlsx_folder}/{filename}', engine='openpyxl', sheet_name='embryo')
        individual_embryo_resolution_data['embryo'] = embryo_name
        embryo_resolution_data.append(individual_embryo_resolution_data)
    return pd.concat(embryo_resolution_data).reset_index(drop=True)

def get_volume_regression_coefficients(embryo_resolution_data):
    """
    # Description
    ---
    Computes the linear regression of an embryo volume and its time points for each unique embryo.
    Returns the coefficients a and b respectively for each embryo in a dict.
    @ Gregoire Malandain

    # Argument(s)
    ---
        `embryo_resolution_data` (pandas df): data structure generated using concatenate_embryo_resolution_data().

    # Usage
    ---
    >>> data = concatenate_embryo_resolution_data('xlsx')
    >>> print(get_volume_regression_coefficients(data))
    >>> print(get_volume_coefficients(data['embryo']))
    ... {'Astec-pm7': (-29805.92158159059, 74225682.50753953), ... 
    'Astec-pm8': (-24884.072590699667, 52644557.01324528)}
    """
    vol_coefficients = {}
    for embryo in embryo_resolution_data['embryo'].unique():
        individual_embryo_data = embryo_resolution_data[embryo_resolution_data['embryo']==embryo]
        xnp = individual_embryo_data['tp'].to_numpy()[:, np.newaxis]
        ynp = individual_embryo_data['volume'].to_numpy()[:, np.newaxis]
        ransac = RANSACRegressor()
        ransac.fit(xnp, ynp)
        vol_coefficients[embryo] = (ransac.estimator_.coef_[0][0], ransac.estimator_.intercept_[0])
    return vol_coefficients

def get_voxelsize_correction(embryo_resolution_data, target_volume=target_volume):
    """
    # Description
    ---
    Returns a dict associating a voxelsize correction to an embryo for each time point.
    The correction is calculated by aligning the embryo volumes to a target value.
    The volume and surface variability can still be explored at the cell- and tissue- resolutions.
    @ Gregoire Malandain.
    
    # Argument(s)
    ---
        `embryo_resolution_data` (pandas df): data structure generated using concatenate_embryo_resolution_data().
        `target_volume` (numerical): value of the embryo volume after alignment.

    # Usage
    ---
    >>> data = concatenate_embryo_resolution_data('xlsx')
    >>> print(get_voxelsize_correction(data)['Astec-pm7'])
    ... {4: 0.23812308016809025, 5: 0.23815449069993436, ... 82: 0.2406240861870154}
    """
    voxelsize_correction = {}
    volume_regression_coefficients = get_volume_regression_coefficients(embryo_resolution_data)
    for embryo in embryo_resolution_data['embryo'].unique():
        voxelsize_correction[embryo] = {}
        individual_embryo_data = embryo_resolution_data[embryo_resolution_data['embryo']==embryo]
        for tp in individual_embryo_data.sort_values(by='tp')['tp']:
            volume_at_t = volume_regression_coefficients[embryo][0] * tp + volume_regression_coefficients[embryo][1]
            voxelsize_correction[embryo][tp] = np.cbrt(target_volume / volume_at_t)
    return voxelsize_correction

def apply_voxelsize_correction(xlsx_folder):
    """
    # Description
    ---
    Applies a voxel size correction for each embryo.
    The volume and surface variability can still be explored at the cell- and tissue- resolutions.
    @ Gregoire Malandain
    
    # Argument(s)
    ---
        `xlsx_folder` (str): name of the folder containing the parsed segmented embryos.
    """
    embryo_resolution_data = concatenate_embryo_resolution_data(xlsx_folder)
    voxelsize_correction = get_voxelsize_correction(embryo_resolution_data)

    for filename in listdir(xlsx_folder):
        embryo_name = filename.split('.')[0]
        individual_data = pd.read_excel(f'{xlsx_folder}/{filename}', engine='openpyxl', sheet_name=None)
    
        for tp, corr in voxelsize_correction[embryo_name].items():
            for sheet in ['embryo', 'tissue', 'cell']:
                individual_data[sheet].volume[individual_data[sheet].tp==tp] *= corr**3
                individual_data[sheet].total_surface[individual_data[sheet].tp==tp] *= corr**2
            for sheet in ['tissue_contacts', 'cell_contacts']:
                individual_data[sheet].surface[individual_data[sheet].tp==tp] *= corr**2
    
        with pd.ExcelWriter(f'{xlsx_folder}/{filename}') as writer:
            for sheet, df in individual_data.items():
                df.to_excel(writer, sheet_name=sheet, index=False)

def find_t(n_cells_embryo_per_time, n_cells):
    """
    # Description
    ---
    @ Gregoire Malandain
    """
    time_points = list(n_cells_embryo_per_time['tp'])

    if n_cells in n_cells_embryo_per_time['cell_count'].values:
        times = [tp for tp in time_points if n_cells_embryo_per_time[n_cells_embryo_per_time.tp==tp].cell_count.values[0]==n_cells]
        return (min(times) + max(times)) / 2.0

    smaller_times = [tp for tp in time_points if n_cells_embryo_per_time.loc[n_cells_embryo_per_time.tp==tp, 'cell_count'].values[0] < n_cells]
    larger_times = [tp for tp in time_points if n_cells_embryo_per_time.loc[n_cells_embryo_per_time.tp==tp, 'cell_count'].values[0] > n_cells]
    
    ts = max(smaller_times)
    tl = min(larger_times)
    n_cells_at_ts = n_cells_embryo_per_time.loc[n_cells_embryo_per_time['tp']==ts, 'cell_count'].values[0]
    n_cells_at_tl = n_cells_embryo_per_time.loc[n_cells_embryo_per_time['tp']==tl, 'cell_count'].values[0]
    return ts + (tl - ts) * (n_cells - n_cells_at_ts) / (n_cells_at_tl - n_cells_at_ts)

def temporal_alignement(reference_cells_per_time, cells_per_time):
    """
    # Description
    ---
    Computes the coefficients required to align an embryo development rate to a reference embryo.
    @ Gregoire Malandain

    # Argument(s)
    ---
        `reference_cells_per_time` (pandas df): a reference embryo dataframe associating a 'cell_count' to a 'tp'.
        `cells_per_time` (pandas df): an embryo dataframe associating a 'cell_count' to a 'tp'.
    """
    ref_cells_count, cells_count = list(reference_cells_per_time.cell_count), list(cells_per_time.cell_count)
    first = max(min(ref_cells_count), min(cells_count))
    last = min(max(ref_cells_count), max(cells_count))
    
    ref_times = []
    times = []
    i = 0
    for n in range(first, last+1):
        if n not in ref_cells_count and n not in cells_count:
            continue
        ref_times += [find_t(reference_cells_per_time, n)]
        times += [find_t(cells_per_time, n)]
        i += 1

    num = sum(np.multiply(times, ref_times)) - sum(times) * sum(ref_times) / len(times)
    den = sum(np.multiply(times, times)) - sum(times) * sum(times) / len(times)
    a = num/den
    b = (sum(ref_times) - a * sum(times)) / len(times)
    
    return a, b 

def get_temporal_coefficients(embryo_resolution_data, reference_embryo_name=reference_embryo['name']):
    """
    # Description
    ---
    Returns a dict associating each unique embryo to two coefficients a and b. 
    These coefficients are used to align the growth rates of individual embryos to a reference.
    @ Gregoire Malandain.
    
    # Argument(s)
    ---
        `embryo_resolution_data` (pandas df): data structure generated using concatenate_embryo_resolution_data().
        `reference_embryo_name` (str): name of the embryo used as reference.

    # Usage
    ---
    >>> data = concatenate_embryo_resolution_data('xlsx')
    >>> print(get_temporal_coefficients(data))
    ... {'Astec-pm7': (0.831013049678967, -9.799705044691862), ... 'Astec-pm8': (1.0, 0.0)}
    """
    coefficients = {}
    reference_cells_per_time = embryo_resolution_data[embryo_resolution_data['embryo']==reference_embryo_name][['tp', 'cell_count']]
    for embryo in embryo_resolution_data.embryo.unique():
        cells_per_time = embryo_resolution_data[embryo_resolution_data['embryo']==embryo][['tp', 'cell_count']]
        coefficients[embryo] = temporal_alignement(reference_cells_per_time, cells_per_time)
    return coefficients

def apply_temporal_alignment(xlsx_folder, reference_embryo=reference_embryo):
    """
    # Description
    ---
    Translates the time points into uniform minutes post-fertilization for each embryo.
    A new column is added on each sheet.
    @ Gregoire Malandain

    # Argument(s)
    ---
        `xlsx_folder` (str): name of the folder containing the parsed segmented embryos.
        `reference_embryo` (dict): metadata about the reference embryo.    
    """
    embryo_resolution_data = concatenate_embryo_resolution_data(xlsx_folder)
    coefficients = get_temporal_coefficients(embryo_resolution_data)

    for filename in listdir(xlsx_folder):
        embryo_name = filename.split('.')[0]
        individual_data = pd.read_excel(f'{xlsx_folder}/{filename}', engine='openpyxl', sheet_name=None)
    
        for sheet in ['embryo', 'tissue', 'tissue_contacts', 'cell', 'cell_contacts']:
            scaled_time_point = (
                individual_data[sheet].tp * coefficients[embryo_name][0] + coefficients[embryo_name][1])
            individual_data[sheet]['minutes_post_fertilization'] = (
                scaled_time_point * reference_embryo['time_point_interval_seconds'] + reference_embryo['start_seconds']) // 60

        with pd.ExcelWriter(f'{xlsx_folder}/{filename}') as writer:
            for sheet, df in individual_data.items():
                df.to_excel(writer, sheet_name=sheet, index=False)