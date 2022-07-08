"""
Main ACoV script called in the command line with positional arguments to accomplish specific tasks.
Command line arguments are described in the main_settings script.
_______________

    + parse: store the anatomical properties of embryos into .xlsx files from .xml segmentation files
    + alignment: applies a voxelsize and a time correction on the embryos
    + matrix: converts the .xlsx files into cell composition matrices
    + 
_______________

@ ASLOUDJ Yanis
07/07/2022
"""


from main_settings import *
from scripts.prog import parse

import argparse


    ### command-line arguments parser ###

argparser = argparse.ArgumentParser(description=argparser_description, formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('task', help=task_help)
argparser.add_argument('-g', '--geometry', action='store_true', help=geometry_help)


    ### ACOB run ###

optional_args = {'parse': ['geometry'], 'align': ['geometry']}
command_line_args = argparser.parse_args()
fun_name = command_line_args.task
fun = globals()[fun_name]
fun_args_names = []
for arg in optional_args[fun_name]:
    fun_args_names.append(arg)
fun_args = []
for arg in fun_args_names:
    fun_args.append(command_line_args.__dict__[arg])
fun(*fun_args)
