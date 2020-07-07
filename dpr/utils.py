from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable


def read_per_line(fname, maxlines=50000):

    array = []
    new_line = []
    nlines = 0

    with open(fname, 'r') as file:
        for line in file:
            # empty line, so we read the whole track and must start a new one
            if line in ['\n', '\r\n']:
                array.append(np.array(new_line, dtype=np.float32))
                new_line = []
            else:
                this_line = line.strip()

                # we have more than one coordinate per line, such as a position xyz
                if len(this_line.split()) > 1:
                    this_line = np.array([float(x) for x in this_line.split()])

                new_line.append(this_line)

            if nlines >= maxlines:
                break
            nlines += 1

    return array


def strip_first_col(fname, delimiter=None):
    with open(fname, 'r') as file:
        for line in file:
            try:
                yield line.split(delimiter, 1)[1]
            except IndexError:
                continue


def strip_header(filename, columns=0):
    header = np.genfromtxt(filename, dtype=str, usecols=columns)
    with open(filename, 'r') as file:
        ncols = len(file.readline().split())
    bundles = np.genfromtxt(filename, usecols=range(1, ncols))
    return bundles, header


# Actual colors I used in the manuscript
blue = np.array([0.6093924, 0.73212757, 0.82249106])
green = np.array([0.22953434, 0.57685998, 0.42976558])


def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def draw_fancy_graph(pval, coords1, coords2, truncated_coords1, truncated_coords2, average1, average2, coord1_label='X', coord2_label='Y',
                     pval_threshold=1., pval_cmap=plt.cm.hot, mean_fiber_cmap=green, bundle_cmap=blue,
                     shadow_cmap='gray', title=None, draw_colorbar=True):

    if title is None:
        title = 'p-values after realignment'

    # This is to draw on a white background
    with plt.style.context('seaborn-whitegrid'):

        fig, ax = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(8, 8))

        # Draw the full length shadow bundle
        for x, y in zip(coords1, coords2):
            ax.plot(x, y, color=shadow_cmap, alpha=0.1, zorder=1)

        # Draw the original coord, but truncated between rois
        for x, y in zip(truncated_coords1, truncated_coords2):
            ax.plot(x, y, color=bundle_cmap, alpha=0.3, zorder=2)

        # Draw the mean coord
        x = average1
        y = average2
        ax.plot(x, y, color=mean_fiber_cmap, zorder=5)

        # We resample the pvals because the coords may be at the tracking stepsize and the final results at the voxel resolution
        old = np.arange(len(pval)) / len(pval)
        new = np.arange(len(x)) / len(x)
        pval_resampled = np.interp(new, old, pval)

        # This makes everything above the threshold invisible on the final graph
        pval_resampled[pval_resampled > pval_threshold] = np.nan
        cmap_axis = ax.scatter(x, y, c=pval_resampled, cmap=pval_cmap, zorder=10, marker='.', vmin=0, vmax=pval_threshold)

        if draw_colorbar:
            colorbar(cmap_axis)

        ax.axis('equal')
        ax.set_xlabel("{} coordinates (mm)".format(coord1_label), fontsize=12)
        ax.set_ylabel("{} coordinates (mm)".format(coord2_label), fontsize=12)

        ax.set_title(title, fontsize=20)
        fig.tight_layout()

    return fig, ax
