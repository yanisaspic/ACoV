"""
Main ACoV script called in the command line with positional arguments to accomplish specific tasks.
Command line arguments are described in the main_settings script.
_______________

    + interactive approach using argparse and the command line
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

    
    ### command-line interpreter ###

def interpret_command_line(argument_parser):
    """
    # Description
    ---
    Parses the arguments on the command line to identify a specific task, and fulfill it.

    # Argument(s)
    ---
        `argument_parser` (argparse.ArgumentParser object).
    """
    command_line_arguments = argument_parser.parse_args()
    task_name = command_line_arguments.task
    task_function = globals()[command_line_arguments.task]

    function_arguments = []
    try:
        for arg in task_specific_optional_arguments[task_name]:
            function_arguments.append(command_line_arguments.__dict__[arg])
    except KeyError:
        pass
    
    task_function(*function_arguments)


    ### ACoV ###

interpret_command_line(argparser)