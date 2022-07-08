"""
Global variables called in the script main.py.
_______________

    + ACoV documentation for the command line help
    + task-specific optional arguments
_______________

@ ASLOUDJ Yanis
07/08/2022
"""


    ### ordered task-specific optional arguments ###

task_specific_optional_arguments = {
    'parse': ['geometry'], 
    'align': ['geometry']}


    ### command-line arguments docstrings ###

argparser_description="""Explore the embryogenesis cell composition variability from segmented embryos.
Using a binary matrix-based approach, multiple figures can be generated to study this variability.
Data pre-processing tasks are available and can be re-used to start studying the geometry variability. :-)
"""
task_help = """{parse, align, matrix, preprocess}
parse       stores the anatomical properties of embryos into multi-resolution .xlsx with human-readable objects and time points.
    *requires .xml segmentation files
preprocess  parse, align and matrix at once.
"""
geometry_help = """{parse, align, preprocess}
[parse] also computes implicit geometrical properties: object total surface, surface and volume ratios.
    *greatly increases computation time
[align] also applies a voxelsize correction of the geometrical properties to allow the comparison between embryos.
"""