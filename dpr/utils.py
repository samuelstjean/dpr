from __future__ import division, print_function

import numpy as np


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
