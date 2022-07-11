"""
Global variables called in the script alignment.py.
_______________

    + target value for the voxelsize correction
    + informations about the reference embryo used for the growth rate alignment
_______________

@ ASLOUDJ Yanis
07/10/2022
"""


    ### voxelsize correction ###

target_volume = 10**6


    ### reference embryo informations ###

# Astec-pm8 segmentation starts 4 hours post-fertilization and a snapshot is taken every 2 minutes :
reference_embryo = {
    'name': 'Astec-pm8',
    'start_seconds': 4 * 3600,
    'time_point_interval_seconds': 2 * 60
    }