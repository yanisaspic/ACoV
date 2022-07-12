"""
API called by the script prog.py to generate figures.
_______________

    + illustrations of the embryo alignments
    + variability dynamics during the embryo development
_______________

@ ASLOUDJ Yanis
07/11/2022
"""


from scripts.alignment import concatenate_unique_resolution_data

import seaborn as sns
import matplotlib.pyplot as plt


def temporal_alignment_figure(xlsx_folder, figure_folder):
    """
    # Description
    ---
    Generates two plots of the embryo cell count, according to the specific time points and the uniform minutes post-fertilization respectively.

    # Argument(s)
    ---
        `xlsx_folder` (str): name of the folder containing the parsed segmented embryos.
        `figure_folder` (str): name of the folder where the figures are saved.
    """
    sns.set()
    fig, axs = plt.subplots(1, 2, figsize=(12,12))
    embryo_resolution_data = concatenate_unique_resolution_data(xlsx_folder)
    time_axes = ['tp', 'minutes_post_fertilization']

    for i in range(len(time_axes)):
        sns.lineplot(data=embryo_resolution_data, x=time_axes[i], y='cell_count', hue='embryo', ax=axs[i])

    fig.suptitle("Temporal alignment of the embryos.", fontsize=24)
    fig.tight_layout()
    plt.savefig(f'{figure_folder}/temporal_alignment.png')

def lineplot_variability(variability_dataframe, x_axis, y_axis, filename, hue=None):
    """
    # Description
    ---
    Generates a lineplot with the variability as y-axis and a proxy of time as the x-axis.
    The input dataframe has time as columns and a modality of variability (e.g. calculation method or tissue subgroup) as rows.
    If the dataframe has multiple rows, each modality is plotted in a unique color.

    # Argument(s)
    ---
        `variability_dataframe` (pandas df): df with time as columns and a modality of variability as row. 
        `filename` (str): path to the generated figure.
    """
    sns.set()
    fig, axs = plt.subplots(figsize=(12,12))
    try:
        sns.lineplot(data=variability_dataframe, x=x_axis, y=y_axis, hue=hue, ax=axs)
    except ValueError:
        sns.lineplot(data=variability_dataframe, x=x_axis, y=y_axis, ax=axs)
    fig.tight_layout()
    plt.savefig(filename)

def plot_heatmap_distances(distances_matrix, filename):
    """
    # Description
    ---
    Plots a heatmap of relative distances.
    """
    fig, axs = plt.subplots(figsize=(12, 12))
    sns.heatmap(distances_matrix, cmap='coolwarm', fmt='g', ax=axs)
    fig.tight_layout()
    plt.savefig(filename)