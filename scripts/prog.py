"""
API called by the script main.py to fulfill the tasks handled by ACoV.
_______________

@ ASLOUDJ Yanis
07/08/2022
"""


from scripts.parser import parse_xml_completely, save_multi_resolution_data
from scripts.alignment import apply_voxelsize_correction, apply_temporal_alignment, concatenate_unique_resolution_data
from scripts.matrix import get_cell_count_matrices, save_cell_count_matrices
from scripts.variability import get_cell_resolution_variability, get_tissue_resolution_variability
from scripts.graphics import *

import pandas as pd
from os import listdir, rename
from time import time


def parse(geometry=False):
    """
    # Description
    ---
    Converts all the .xml files stored in an input folder into .xlsx files. They are stored in an output folder. 
    Afterwards, the input files are also moved to a throwaway folder.

    # Argument(s)
    ---
        `geometry` (bool): if True, updates the geometrical properties of surface and volume (time-consuming).
    """
    start = time()
    xml_files = listdir('data/xml')
    total = len(xml_files)
    print(f'{total} .xml files ready to be parsed:')

    for filename in xml_files:
        file_start = time()
        print(f'\tCurrently parsing {filename}...')
        embryo_name = filename.split('.')[0]
        data = parse_xml_completely(f'data/xml/{filename}', geometry)
        save_multi_resolution_data(data, f'data/xlsx/{embryo_name}.xlsx')
        rename(f'data/xml/{filename}', f'data/used/{filename}')
        print(f'\t\t-> done in {round(time() - file_start)} seconds.')
    
    print(f'{total} .xml files parsed in {round(time() - start)} seconds.')

def align(geometry=False):
    """
    # Description
    ---
    Updates the .xlsx files to facilitate the comparison between different embryos.

    # Argument(s)
    ---
        `geometry` (bool): if True, applies a voxelsize correction to the geometrical properties (time-consuming).
    """
    start = time()
    xlsx_files = listdir('data/xlsx')
    total = len(xlsx_files)
    print(f'{total} .xlsx files ready to be aligned:')

    if geometry:
        start_voxelsize = time()
        print('\tAlignment of the voxelsize...')
        apply_voxelsize_correction('data/xlsx')
        print(f'\t\t-> done in {round(time() - start_voxelsize)} seconds.')
    
    start_temporal = time()
    print('\tAlignment of the growth rate...')
    apply_temporal_alignment('data/xlsx')
    print(f'\t\t-> done in {round(time() - start_temporal)} seconds.')

    print(f'{total} .xlsx files aligned in {round(time() - start)} seconds.')
    temporal_alignment_figure('data/xlsx', 'figures')

def matrix():
    """
    # Description
    ---
    Generates cell composition matrices from the .xlsx files.
    """
    start = time()
    xlsx_files = listdir('data/xlsx')
    total = len(xlsx_files)
    print(f'{total} .xlsx files ready to be matrixed:')

    for resolution in ['cell', 'tissue']:
        sub_resolution_data = concatenate_unique_resolution_data('data/xlsx', resolution)
        matrices = get_cell_count_matrices(sub_resolution_data)
        for label, time_axis in {'ecc': 'embryo_cell_count', 'mpf': 'minutes_post_fertilization'}.items():
            save_cell_count_matrices(matrices[time_axis], f'data/mat/{resolution}_{label}')

    print(f'{total} .xlsx files matrixed in {round(time() - start)} seconds.')

def preprocess(geometry=False):
    """
    # Description
    ---
    Completes the parse, align and matrix tasks consecutively.

    # Argument(s)
    ---
        `geometry` (bool): cf. parse or align.
    """
    parse(geometry)
    align(geometry)
    matrix()

def plot_cell_resolution_variability(cell_matrices_dict, time_axis='mpf'):
    """
    # Description
    ---
    Generates plots of variability from cell composition matrices at the cell-resolution :
        + lineplots of the min, max and mean number (or percentage) of unique cells during the embryo development
        + lineplot of the evenness indices dynamics during the embryo development
        + heatmap of the variability based on the pairwise distance computation method

    # Argument(s)
    ---
        `cell_matrices_dict` (dict): embryo names as keys and cell composition matrices at the cell-resolution as values.
        `time_axis` (str): 'ecc' or 'mpf'.
    """
    embryo_cell_count, variability_axis = False, 'cityblock distance'
    if time_axis == 'ecc':
        embryo_cell_count, variability_axis = True, 'percentage of cells unique to an embryo'
    cell_resolution_variability = get_cell_resolution_variability(cell_matrices_dict, embryo_cell_count)

    cell_resolution_variability['evenness'] = cell_resolution_variability['evenness'].melt(ignore_index=False).reset_index()
    cell_resolution_variability['evenness'].columns = ['method', time_axis, 'evenness']
    lineplot_variability(cell_resolution_variability['evenness'], time_axis, 'evenness', f'figures/cell/{time_axis}_evenness.png', hue='method')

    for modality in ['min', 'max', 'mean']:
        cell_resolution_variability[modality] = cell_resolution_variability[modality].melt()
        cell_resolution_variability[modality].columns = [time_axis, variability_axis]
        lineplot_variability(cell_resolution_variability[modality], time_axis, variability_axis, f'figures/cell/{time_axis}_{modality}.png')

    plot_heatmap_distances(cell_resolution_variability['distances'], f'figures/cell/{time_axis}_distances.png')

def plot_tissue_resolution_variability(tissue_matrices_dict, time_axis='mpf'):
    """
    # Description
    ---
    Generates plots of variability from cell composition matrices at the tissue-resolution :
        + lineplots of the min, max and mean number (or percentage) of cells uniquely located in a tissue during the embryo development
        + similar lineplot for specific groups of tissues

    # Argument(s)
    ---
        `cell_matrices_dict` (dict): embryo names as keys and cell composition matrices at the cell-resolution as values.
        `time_axis` (str): 'ecc' or 'mpf'.
    """
    
        ### Python weird behavior (cf. get_tissue_resolution_variability)

    # if time_axis == 'ecc':
    #     embryo_cell_count, variability_axis = True, 'percentage of cells uniquely located in a tissue of an embryo'

    embryo_cell_count, variability_axis = False, 'cityblock distance'
    tissue_resolution_variability = get_tissue_resolution_variability(tissue_matrices_dict, embryo_cell_count)

    tissue_resolution_variability['subgroup'] = tissue_resolution_variability['subgroup'].melt(ignore_index=False).reset_index()
    tissue_resolution_variability['subgroup'].columns = ['tissues', time_axis, variability_axis]
    lineplot_variability(tissue_resolution_variability['subgroup'], time_axis, variability_axis, f'figures/tissue/{time_axis}_subgroup.png', hue='tissues')

    for modality in ['min', 'max', 'mean']:
        tissue_resolution_variability[modality] = tissue_resolution_variability[modality].melt()
        tissue_resolution_variability[modality].columns = [time_axis, variability_axis]
        lineplot_variability(tissue_resolution_variability[modality], time_axis, variability_axis, f'figures/tissue/{time_axis}_{modality}.png')
    
def plot_variability(matrices_filename):
    """
    # Description
    ---
    Generates plots of variability from cell composition matrices. The resulting plots are dependent on the resolution.

    # Argument(s)
    ---
        `matrices_filename` (str): name of the .xlsx matrices file.
    """
    callback = {'cell': plot_cell_resolution_variability, 'tissue': plot_tissue_resolution_variability}

    data = pd.read_excel(matrices_filename, engine='openpyxl', sheet_name=None, index_col=0)
    [resolution, time_axis] = matrices_filename.split('.')[0].split('/')[-1].split('_')

    return callback[resolution](data, time_axis)

def plot():
    """
    # Description
    ---
    From the cell composition matrices geneared with matrix(), plot figures to explore the division variability.
    """
    start = time()
    print(f'Starts to plot...')
    for file in listdir('data/mat'):
        plot_variability(f'data/mat/{file}')
    print(f'Plotting done in {round(time() - start)} seconds.')