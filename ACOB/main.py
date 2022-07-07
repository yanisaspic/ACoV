"""
Main script of ACOB called in the command line with a keyword argument :
_______________

    + parse: converts the .xml embryo files into multi-resolution .xlsx files
    + alignment: applies a voxelsize and a time correction on the embryos
    + matrix: converts the .xlsx files into cell composition matrices
    + 
_______________

@ ASLOUDJ Yanis
07/07/2022
"""

from scripts.parser import parse_xml_completely, save_multi_resolution_data

import argparse
from os import listdir, rename


def parse():
    """
    # Description
    ---
    Parses all the .xml files found inside a directory and stores the data .
    """
    for filename in listdir('data/XML'):
        embryo_name = filename.split('.')[0]
        data = parse_xml_completely(f'data/XML/{filename}')
        save_multi_resolution_data(data, f'data/XLSX/{embryo_name}.xlsx')
        rename(f'data/XML/{filename}', f'data/XML_parsed/{filename}')


argument_parser = argparse.ArgumentParser(description="""Explore the embryogenesis cell composition variability from segmented embryos.""")
argument_parser.add_argument('-task', help='indicate a task to fulfill: parse, alignment, matrix')
args = argument_parser.parse_args()
fun = globals()[args.task]
fun()