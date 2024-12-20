import numpy as np

import argparse
import logging
import os

from dpr.utils import draw_fancy_graph, read_per_line
from ast import literal_eval


DESCRIPTION = """
Script to draw overlays after applying the diffusion profile realignment (DPR) algorithm.
The x and y axis relate to the graph itself, but any spatial combination can be used.
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

    p.add_argument('pvals',
                   help='Path of the input text file to overlay.')

    p.add_argument('coords',
                   help='Path of the input text file for the whole coordinates.')

    p.add_argument('truncated',
                   help='Path of the input text file for the x axis coordinates between rois.')

    p.add_argument('representative',
                   help='Path of the input text file providing the representative streamline on the x axis.')

    p.add_argument('columns', type=literal_eval,
                   help='The two columns (0, 1 or 2) to use when reading from the text files containing xyz coordinates.\n'
                   'Separate them by a comma (without space) like this : 0,1')

    p.add_argument('output', metavar='output',
                   help='Path of the output figure.')

    p.add_argument('--labelx', metavar='string', default='X',
                   help='Label for the x axis.')

    p.add_argument('--labely', metavar='string', default='Y',
                   help='Label for the y axis.')

    p.add_argument('--title', metavar='string',
                   help='Title of the figure.')

    p.add_argument('--dpi', metavar='int', default=150, type=int,
                   help='Dot per inch (resolution) of the figure (default : 150)')

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

    axis1, axis2 = args.columns

    # Load up the simple text files
    pvals = np.loadtxt(args.pvals)

    representative = np.loadtxt(args.representative)
    representative1 = representative[:, axis1]
    representative2 = representative[:, axis2]

    # Load up the text files which are lists of various points
    coords = read_per_line(args.coords)
    coords1 = [coords[i][:, axis1] for i in range(len(coords))]
    coords2 = [coords[i][:, axis2] for i in range(len(coords))]

    truncated = read_per_line(args.truncated)
    truncated1 = [truncated[i][:, axis1] for i in range(len(truncated))]
    truncated2 = [truncated[i][:, axis2] for i in range(len(truncated))]

    # Finally draw up everything
    fig, axes = draw_fancy_graph(pvals, coords1, coords2, truncated1, truncated2, representative1, representative2,
                                 coord1_label=args.labelx, coord2_label=args.labely, title=args.title)

    fig.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
