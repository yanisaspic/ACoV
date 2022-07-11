"""
API called by the script main.py to fulfill the tasks handled by ACoV.
_______________

@ ASLOUDJ Yanis
07/08/2022
"""


from scripts.parser import parse_xml_completely, save_multi_resolution_data
from scripts.alignment import apply_voxelsize_correction, apply_temporal_alignment, concatenate_unique_resolution_data
from scripts.matrix import get_cell_count_matrices, save_cell_count_matrices
from scripts.graphics import *

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
        print(f'\t\t-> done in {round(time() - start_voxelsize())} seconds.')
    
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