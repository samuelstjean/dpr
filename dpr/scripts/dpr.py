import numpy as np
import matplotlib.pyplot as plt

import argparse
import logging
import os

from dpr.register import align_bundles, resample_bundles_to_same, flip_fibers, truncate

# These are a few reading functions for text file I made, feel free to use them if they fit your data
from dpr.utils import read_per_line, strip_first_col, strip_header


DESCRIPTION = """
Main script for the diffusion profile realignment (DPR) algorithm.
"""

EPILOG = """
Reference :

St-Jean, Samuel, Chamberland, Maxime, Viergever, Max A. and Leemans, Alexander.
Reducing variability in along-tract analysis with diffusion profile realignment, NeuroImage, 2019, ISSN 1053-8119

Available at: https://doi.org/10.1016/j.neuroimage.2019.06.016
"""


def buildArgsParser():

    p = argparse.ArgumentParser(description=DESCRIPTION,
                                epilog=EPILOG,
                                formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument('data', metavar='input',
                   help='Path of the input text file of bundles.')

    p.add_argument('output', metavar='output',
                   help='Path of the output text file of realigned bundles.')

    p.add_argument('--exploredti', action='store_true',
                   help='Strip the first column from the input file as used by explore dti to store each subject name.')

    p.add_argument('--do_graph', action='store_true',
                   help='Save a small plot of the original and realigned data.')

    p.add_argument('--points', metavar='int', type=int,
                   help='Number of points for the final resampling. Default: longest bundle')

    p.add_argument('--coordinates', metavar='input',
                   help='Text file with the xyz coordinates of each bundle to determine the initial system of coordinates.'
                   '\nUseful if your data is not already in increasing order.')

    p.add_argument('-f', '--force', action='store_true', dest='overwrite',
                   help='If set, overwrites the output text file if it already exists.')

    p.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                   help='If set, print useful information during processing.')

    p.add_argument('-l', '--log', dest='logfile', metavar='file',
                   help='Save the logging output to this file. Implies verbose output.')

    return p


def main():
    parser = buildArgsParser()
    args = parser.parse_args()

    logger = logging.getLogger('diffusion profile realignment')

    if args.logfile is not None:
        handler = logging.FileHandler(args.logfile)
        args.verbose = True
    else:
        handler = logging.StreamHandler(args.logfile)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M:%S")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if args.verbose:
        logger.setLevel(logging.INFO)
        logger.info('Verbosity is on')

    if os.path.isfile(args.output):
        if args.overwrite:
            logger.warning(f'Overwriting {os.path.realpath(args.output)}')
        else:
            parser.error(f'{args.output} already exists! Use -f or --force to overwrite it.')

    # load the data
    _, ext = os.path.splitext(args.data)

    if ext == '.csv':
        delimiter = ','
    else:
        delimiter = None

    if args.exploredti:
        data, header = strip_header(args.data, delimiter=delimiter)
    else:
        data = np.genfromtxt(args.data, delimiter=delimiter)

    logger.info(f'Number of bundles is {data.shape[0]}, original number of points per tract is {data.shape[1]}')

    if args.coordinates is not None:
        coordinates = np.genfromtxt(args.coordinates)
        data = flip_fibers(data, coordinates)

    aligned, shifts = align_bundles(data)
    truncated = truncate(aligned)

    # resample so that the number of points is the same as the original file or as requested
    if args.points is None:
        num_points = truncated.shape[1]
    else:
        num_points = args.points

    resampled = resample_bundles_to_same(truncated, num_points=num_points)

    logger.info(f'Final number of points per tract is {num_points}')

    if args.do_graph:
        f, axes = plt.subplots(2, 1, sharex=False, sharey=True)
        labels = 'Before realignment', 'After realignment'

        for bund, ax, label in zip([data, resampled], axes, labels):
            mean = np.nanmean(bund, axis=0)
            std = np.nanstd(bund, axis=0)
            ax.fill_between(range(len(mean)), mean - std, mean + std)
            ax.plot(mean, color='r', label=label)

            ax.set_xlim(0, len(mean))
            ax.set_ylim(0, None)
            ax.legend(loc='lower left', fontsize=12)
            ax.set(ylabel='Metric', xlabel='Coordinates')

        f.tight_layout()
        f.savefig(args.output.replace('.txt','.png'), dpi=100, bbox_inches='tight')

    # should patch back whatever was removed beforehand here
    if args.exploredti:
        resampled = np.column_stack((header, resampled))

    _, ext = os.path.splitext(args.output)

    if ext == '.csv':
        delimiter = ','
    else:
        delimiter = ' '

    with open(args.output, 'w') as file:
        np.savetxt(file, resampled, delimiter=delimiter, fmt='%s')
