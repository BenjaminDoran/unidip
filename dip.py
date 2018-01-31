"""
    Module for computing The Hartigans' dip statistic

    The dip statistic measures unimodality of a sample from a random process.

    See:
    Hartigan, J. A.; Hartigan, P. M. The Dip Test of Unimodality. The Annals
    of Statistics 13 (1985), no. 1, 70--84. doi:10.1214/aos/1176346577.
    http://projecteuclid.org/euclid.aos/1176346577.
"""

import collections
import numpy as np

def _gcm_(cdf, idxs):
    work_cdf = cdf
    work_idxs = idxs
    gcm = [work_cdf[0]]
    touchpoints = [0]
    while len(work_cdf) > 1:
        distances = work_idxs[1:] - work_idxs[0]
        slopes = (work_cdf[1:] - work_cdf[0]) / distances
        minslope = slopes.min()
        minslope_idx = np.where(slopes == minslope)[0][0] + 1
        gcm.extend(work_cdf[0] + distances[:minslope_idx] * minslope)
        touchpoints.append(touchpoints[-1] + minslope_idx)
        work_cdf = work_cdf[minslope_idx:]
        work_idxs = work_idxs[minslope_idx:]
    return np.array(np.array(gcm)), np.array(touchpoints)

def _lcm_(cdf, idxs):
    g, t = _gcm_(1-cdf[::-1], idxs.max() - idxs[::-1])
    return 1-g[::-1], len(cdf) - 1 - t[::-1]

def _touch_diffs_(part1, part2, touchpoints):
    diff = np.abs((part2[touchpoints] - part1[touchpoints]))
    return diff.max(), diff

def diptst(idxs, nboot):
    """ diptest with pval """
    # sample dip
    d, (_, idxs, left, _, right, _) = dip(idxs=idxs)

    # simulate from null uniform
    unifs = np.random.uniform(size=nboot * idxs.shape[0])\
                     .reshape([nboot, idxs.shape[0]])
    unif_dips = np.apply_along_axis(dip, 1, unifs, just_dip=True)


    # count dips greater or equal to d, add 1/1 to prevent a pvalue of 0
    p = None if unif_dips.sum() == 0\
        else (np.less(d, unif_dips).sum() + 1) / (np.float(nboot) + 1)

    return (d, p, (idxs[len(left)], idxs[-len(right)]))

def dip(idxs=None, histogram=None, just_dip=False):
    """
        Compute the Hartigans' dip statistic either for a histogram of
        samples (with equidistant bins) or for a set of samples.
    """
    if idxs is None:
        idxs = np.arange(len(histogram))
    elif histogram is None:
        h = collections.Counter(idxs)
        idxs = np.msort(list(h.keys()))
        histogram = np.array([h[i] for i in idxs])
    else:
        if len(histogram) != len(idxs):
            raise ValueError("Need exactly as many indices as histogram bins.")
        if len(idxs) != len(set(idxs)):
            raise ValueError("idxs must be unique if histogram is given.")
        if not np.array_equal(np.msort(idxs), idxs):
            idxs_s = np.argsort(idxs)
            idxs = np.asarray(idxs)[idxs_s]
            histogram = np.asarray(histogram)[idxs_s]

    cdf = np.cumsum(histogram, dtype=float)
    cdf /= cdf[-1]

    # check for case 1<N<4 or all identical values
    if len(idxs) <= 4 or idxs[0] == idxs[-1]:
        left = []
        right = [1]
        d = 0.0
        return d if just_dip else (d, (cdf, idxs, left, None, right, None))

    work_idxs = idxs
    work_histogram = np.asarray(histogram, dtype=float) / np.sum(histogram)
    work_cdf = cdf

    D = 0
    left = [0]
    right = [1]

    while True:
        left_part, left_touchpoints = _gcm_(work_cdf-work_histogram, work_idxs)
        right_part, right_touchpoints = _lcm_(work_cdf, work_idxs)

        d_left, left_diffs = _touch_diffs_(left_part,
                                           right_part, left_touchpoints)
        d_right, right_diffs = _touch_diffs_(left_part,
                                             right_part, right_touchpoints)

        if d_right > d_left:
            xr = right_touchpoints[d_right == right_diffs][-1]
            xl = left_touchpoints[left_touchpoints <= xr][-1]
            d = d_right
        else:
            xl = left_touchpoints[d_left == left_diffs][0]
            xr = right_touchpoints[right_touchpoints >= xl][0]
            d = d_left

        left_diff = np.abs(left_part[:xl+1] - work_cdf[:xl+1]).max()
        right_diff = np.abs(right_part[xr:]
                            - work_cdf[xr:]
                            + work_histogram[xr:]).max()

        if d <= D or xr == 0 or xl == len(work_cdf):
            the_dip = max(np.abs(cdf[:len(left)] - left).max(),
                          np.abs(cdf[-len(right)-1:-1] - right).max())
            if just_dip:
                return the_dip/2
            else:
                return the_dip/2, (cdf, idxs, left, left_part, right, right_part)
        else:
            D = max(D, left_diff, right_diff)

        work_cdf = work_cdf[xl:xr+1]
        work_idxs = work_idxs[xl:xr+1]
        work_histogram = work_histogram[xl:xr+1]

        left[len(left):] = left_part[1:xl+1]
        right[:0] = right_part[xr:-1]

def plot_dip(idxs=None, histogram=None):
    """ Plot ECDF with colored intervals """
    from matplotlib import pyplot as plt

    d, (cdf, idxs, left, left_part, right, right_part) = dip(histogram, idxs)


    plt.plot(idxs[:len(left)], left, color='red')
    plt.plot(idxs[len(left)-1:len(left)+len(left_part) - 1], left_part, color='gray')
    plt.plot(idxs[-len(right):], right, color='blue')
    plt.plot(idxs[len(cdf)-len(right)+1-len(right_part):len(cdf)-len(right)+1],
             right_part, color='gray')

    the_dip = max(np.abs(cdf[:len(left)] - left).max(),
                  np.abs(cdf[-len(right)-1:-1] - right).max())
    l_dip_idxs = np.abs(cdf[:len(left)] - left) == the_dip
    r_dip_idxs = np.abs(cdf[-len(right)-1:-1] - right) == the_dip
    print(the_dip/2, d)

    plt.vlines(x=idxs[:len(left)][l_dip_idxs],
               ymin=cdf[:len(left)][l_dip_idxs],
               ymax=cdf[:len(left)][l_dip_idxs] - the_dip)
    plt.vlines(x=idxs[-len(right):][r_dip_idxs],
               ymin=cdf[-len(right):][r_dip_idxs],
               ymax=cdf[-len(right):][r_dip_idxs] + the_dip)

    plt.plot(np.repeat(idxs, 2)[1:], np.repeat(cdf, 2)[:-1], color='black')
    plt.scatter(idxs, cdf)

    plt.show()


def crit_points(random_function, quantiles, sample_size, n_samples):
    """
        Compute the quantiles for the dip statistic for n_samples
        samples of size sample_size from the random process given by
        random_function.

        Parameters:
        random_function : a paramter-free function which returns random values.
        quantiles : a sequence of values between 0 and 1
        sample_size : the size of the samples to draw from random_function
        n_samples : the number of samples for which to compute dips

        Returns: a list such that the i'th value is the greatest dip observed
        such that the fraction of dips less than or equal to that value is less
        than the i'th value from quantiles.
    """
    data = [[random_function() for _ in range(sample_size)] for __ in range(n_samples)]
    dips = np.array([dip(idxs=samples)[0] for samples in data])

    return np.percentile(dips, [p * 100 for p in quantiles])
