from __future__ import division, print_function

import numpy as np

from dpr.register import truncate, extrapolate


def test_extrapolate():

    x = np.arange(-10, 10)
    a, b, c = -1, 0, 4
    y = a * x**2 + b * x + c
    ymax = y.argmax()

    peak, value = extrapolate([x[ymax - 1], x[ymax], x[ymax + 1]],
                              [y[ymax - 1], y[ymax], y[ymax + 1]],
                              return_value=True)

    np.testing.assert_array_almost_equal(peak, 0)
    np.testing.assert_array_almost_equal(value, 4)


def test_truncate():
    a = np.full((1, 10), np.nan)
    b = np.full((1, 10), np.nan)
    a[:, 2:6] = range(2, 6)
    b[:, 4:7] = range(4, 7)

    c = np.concatenate((a, b))

    answer = np.array([[4, 5], [4, 5]])
    np.testing.assert_equal(truncate(c, mode='shortest'), answer)

    answer = np.zeros((2, len(range(2, 7)))) * np.nan
    answer[0, range(4)] = range(2, 6)
    answer[1, range(2, 5)] = range(4, 7)
    np.testing.assert_equal(truncate(c, mode='longest'), answer)
