"""
API called by the script main.py to fulfill the tasks handled by ACoV.
_______________

    + parse: store the anatomical properties of embryos into .xlsx files from .xml segmentation files
_______________

@ ASLOUDJ Yanis
07/08/2022
"""


from scripts.parser import parse_xml_completely, save_multi_resolution_data

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