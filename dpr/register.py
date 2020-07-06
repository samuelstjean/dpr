from __future__ import division, print_function

import numpy as np

from itertools import product
from warnings import warn

from scipy import ndimage


def align_bundles(bundles, percent=15, padding=0., order=1, eps=1e-5, mode='full_template',
                  remove_outliers=True, remove_baseline=True, whiten=True, normalize=False,
                  return_shifts_matrix=False, rematch_outliers=True):

    bundles = np.array(bundles, copy=True)

    # if we zero padded, put them to nan for now
    bundles[bundles == padding] = np.nan

    npairs = len(bundles)
    size_overlaps = np.isfinite(bundles).sum(axis=1)
    maxoverlaps = np.ceil(percent * size_overlaps / 100)
    shifts = np.zeros((npairs, npairs))

    pairs = filter_pairs(npairs, mode=mode)
    ffts = get_ffts(bundles, whiten=whiten, remove_baseline=remove_baseline)
    ffta = np.zeros_like(ffts[0])
    fftb = np.zeros_like(ffts[0])

    for (a, b) in pairs:
        ffta[:] = ffts[a]
        fftb[:] = ffts[b]
        shifts[a, b] = get_shift_from_fft(ffta, fftb, normalize=normalize)

    # if we used a predefined template, we can just exit now with it
    if len(pairs) < bundles.shape[0]:

        template = mode
        ys = apply_shift(bundles, shifts[template])

        if return_shifts_matrix:
            return ys, shifts, template
        else:
            return ys, shifts[template]

    shifts[np.abs(shifts) < eps] = 0.  # force symmetry at ~zero shift

    # this is how many bundles are inside the overlapping threshold for bundle i
    sumstuff = np.sum(np.abs(shifts) < maxoverlaps, axis=1)

    # if we have more than one template, we pick the one with the min maximum shift
    # this does not seem to happen on real data, only on synthetic
    original_templater = np.argwhere(sumstuff == np.max(sumstuff))

    if len(original_templater) > 1:
        # abs(shifts) will be a symmetric matrix)
        largest_shifts = np.max(np.abs(shifts), axis=1)[original_templater]
        original_templater = original_templater[largest_shifts.argmin()]

    original_templater = int(original_templater)  # it has a weird shape now if singleton

    if rematch_outliers:
        condition = True
        current_outliers = [None]

        # while loop ends when outliers will be empty
        while condition:

            # this is the indices of the outliers for the best alignment template bundle we just found
            outliers = np.abs(shifts[original_templater]) > maxoverlaps
            outliers = np.arange(npairs)[outliers]

            candidates = np.abs(shifts[original_templater]) <= maxoverlaps
            candidates = np.arange(npairs)[candidates]

            for outlier in outliers:

                # find the match that minimize distance to original template without having very large displacement
                new_shifts = shifts[original_templater, :] + shifts[:, outlier]

                # remove candidate templates which are also outliers themselves
                indexes = [idx for idx in np.abs(new_shifts).argsort() if idx in candidates]

                if np.any(indexes):
                    new_shifter = indexes[0]
                    # this is the shift required to realign bundle -> new_template -> original template
                    shifts[outlier, original_templater] = shifts[outlier, new_shifter] + shifts[new_shifter, original_templater]
                    shifts[original_templater, outlier] = shifts[original_templater, new_shifter] + shifts[new_shifter, outlier]

            # We have some outliers which do not overlap between the threshold
            # with the others, so we must break out of an infinite loop.
            # We use a set since order is not important and we replace them outside the loop.
            if set(current_outliers) == set(outliers):
                rematch_outliers = False
                break

            current_outliers = outliers
            condition = len(outliers)
    else:
        # no rematch outliers? put them to zero/nan then
        outliers = np.abs(shifts[original_templater]) > maxoverlaps
        outliers = np.arange(npairs)[outliers]

        if remove_outliers:
            shifts[outliers] = np.nan
            shifts[:, outliers] = np.nan
        else:
            shifts[outliers] = 0.
            shifts[:, outliers] = 0.

        warn('Possible outliers found {}'.format(outliers))

    # use the custom shift patterns
    ys = apply_shift(bundles, shifts[original_templater])

    if return_shifts_matrix:
        return ys, shifts, original_templater
    else:
        return ys, shifts[original_templater]


def resample_bundles_to_same(bundles, num_points=None):

    if bundles.ndim == 1:
        bundles = bundles[None, :]

    if bundles.ndim != 2:
        error = 'bundles needs to be a 2D array, but is {}D'.format(bundles.ndim)
        raise ValueError(error)

    if num_points is None:
        num_points = bundles.shape[1]

    resampled = np.zeros((bundles.shape[0], num_points))

    for i in range(len(bundles)):
        strip = bundles[i]
        length = len(strip)
        old_coord = np.linspace(1, length + 1, num=length, endpoint=True) / length
        new_coord = np.linspace(1, length + 1, num=num_points, endpoint=True) / length
        interp = np.interp(new_coord, old_coord, strip)
        resampled[i] = interp

    return resampled


def apply_shift(bundles, shifts, order=1, padding=np.nan):

    shifted_bundles = np.zeros((bundles.shape[0], 3 * bundles.shape[1]), dtype=np.float32)

    for idx in range(bundles.shape[0]):
        bundle = bundles[idx]
        shift = shifts[idx]

        # if a shift is nan, it's an outlier anyway
        if np.isnan(shift):
            shift = 0

        invoxel_shift, integer_shift = np.modf(shift)
        integer_shift = int(integer_shift)
        pad = np.full(len(bundle), padding)

        new = np.concatenate((pad, bundle, pad))
        new = np.roll(new, integer_shift)
        new = ndimage.shift(new, invoxel_shift, order=order, mode='constant', cval=padding)

        shifted_bundles[idx] = new

    return shifted_bundles


def flip_fibers(bundles, coordinates, padding=np.nan, template=None):
    '''
    bundle - 2D array of M bundles with each metrics of size N as a column
    coordinates - list of points of size whatever x 3
    template - Use this streamline to set the coordinate system.
        If not set, we use the first one from coordinates.
    '''

    if coordinates[0].shape[1] != 3:
        raise ValueError('coordinates should be Nx3 but is {}'.format(coordinates[0].shape))

    new_bundles = np.zeros_like(bundles)

    if template is None:
        template = coordinates[0]

    for idx in range(bundles.shape[0]):
        A = np.sqrt(np.sum((template[0]  - coordinates[idx][0])**2)) + np.sqrt(np.sum((template[-1] - coordinates[idx][-1])**2))
        B = np.sqrt(np.sum((template[-1] - coordinates[idx][0])**2)) + np.sqrt(np.sum((template[-1] - coordinates[idx][0])**2))

        if A > B:
            new_bundles[idx] = np.nan_to_num(bundles[idx][::-1])
            # we still want the padding at the right end though
            cut = np.trim_zeros(new_bundles[idx], trim='f')
            new_bundles[idx] = np.concatenate((cut, np.full(bundles.shape[1] - len(cut), padding)))
        else:
            new_bundles[idx] = bundles[idx]

    return new_bundles


def truncate(bundles, mode='shortest', trimval=np.nan, axis=0):

    if bundles.ndim > 2:
        error = 'Number of dimension must be lower than 2, but was {}'.format(bundles.ndim)
        raise ValueError(error)

    if mode == 'shortest':
        threshold = bundles.shape[0]
    elif mode == 'longest':
        threshold = 1
    elif isinstance(mode, int):

        if mode > 100 or mode < 1:
            raise ValueError('mode must be between 1 and 100 {}'.format(mode))

        threshold = np.floor(mode * bundles.shape[0] / 100)
    else:
        raise ValueError('Unrecognized truncation mode {}'.format(mode))

    if np.isnan(trimval):
        indexes = np.isfinite(bundles).sum(axis=axis) >= threshold
    else:
        indexes = np.sum(bundles != trimval, axis=axis) >= threshold

    truncated = bundles[:, indexes].copy()

    return truncated


def filter_pairs(allpairs, mode):

    if isinstance(allpairs, np.integer) or isinstance(allpairs, int):
        allpairs = list(range(allpairs))

    if isinstance(mode, np.integer) or isinstance(mode, int):
        return [x for x in product([mode], allpairs) if x[1] != mode]
    elif mode == 'full_template':
        pairs = []
        for template in allpairs:
            pairs += [x for x in product([template], allpairs)]
        return pairs
    elif mode == 'half':
        func = np.greater
    elif mode == 'full':
        func = np.not_equal
    elif mode == 'everything':
        func = lambda *_: True  # this always return True
    else:
        error = 'String mismatch, mode was {} of type {}'.format(mode, type(mode))
        raise ValueError(error)

    output = []
    for i, j in product(allpairs, repeat=2):
        if func(j, i):
            output += [(i, j)]

    return output


def get_shift_from_fft(x, y, normalize=False):

    spectrum = crosscorr(x, y, normalize=normalize)

    peak = spectrum.argmax()
    max_value = spectrum[peak]
    # this works because we pad up to 2**N, so it's always even.
    half_length = len(spectrum) // 2

    m0 = spectrum[peak - 1]
    m1 = spectrum[peak + 1]
    p0 = peak - 1
    p1 = peak + 1

    extrapol = extrapolate([p0, peak, p1], [m0, max_value, m1])
    shift = extrapol - half_length
    return shift


def get_ffts(bundles, whiten=True, remove_baseline=True):

    bundles = np.copy(bundles)
    finite = np.isfinite(bundles)

    if remove_baseline:
        for i in range(bundles.shape[0]):
            bundle = bundles[i, finite[i]]
            line = np.linspace(0., np.sqrt(bundle.max()), len(bundle))
            lin_fit = np.polyfit(line, bundle, 1)
            trend = np.polyval(lin_fit, line)
            bundles[i, finite[i]] -= trend

    if whiten:
        for i in range(bundles.shape[0]):
            mean = np.mean(bundles[i, finite[i]])
            std = np.std(bundles[i, finite[i]])
            bundles[i, finite[i]] = (bundles[i, finite[i]] - mean) / std

    N = 2 * bundles.shape[1] - 1
    pad = 2**int(np.ceil(np.log2(N)))
    ffts = np.zeros((bundles.shape[0], pad//2 + 1), dtype=np.complex128)

    for i in range(bundles.shape[0]):
        bundle = bundles[i, finite[i]]
        ffts[i] = np.fft.rfft(bundle, n=pad)

    return ffts


def crosscorr(ffta, fftb, normalize=False):

    A = ffta
    B = fftb.conj()

    if normalize:
        norm = np.abs(A * B)
    else:
        norm = 1

    result = np.fft.irfft(A * B / norm)

    return np.fft.fftshift(result)


def extrapolate(x, y, return_value=False):

    parabola = np.polyfit(x, y, 2)
    a, b, c = parabola

    if a == 0 and b == 0 and c == 0:
        peak = np.nan
        value = np.nan
    else:
        peak = -b / (2 * a)
        value = a * peak**2 + b * peak + c

    if return_value:
        return peak, value
    return peak
