#!python
"""
Author: Benjamin Doran
Date: Jan 2018

Algorithm ref: ...
"""
import numpy as np
import dip

NTRIALS = 100


def unidip(dat, alpha=0.05, _is_model=True, plotdat=None):
    """ Perform unidip algorithm on 1d array """
    intervals = set()
    xs = np.sort(dat)
    _, p, modint = dip.diptst(xs, NTRIALS)

    if not plotdat is None: # if plotting -> show intervals
        _db_plot(plotdat, (xs[0], xs[-1]), intervals | set([modint]))

    # is unimodal, return interval
    if p is None or p > alpha:
        if _is_model:
            intervals.add((xs[0], xs[-1]))
        else:
            wide = _get_full_interval(xs, modint)
            intervals.add(wide)
        return intervals

    # recurse into model interval
    modidx = get_indices(xs, modint)
    intervals |= unidip(xs[modidx[0]:modidx[1]+1], alpha, True, plotdat=plotdat)

    # upper and lower bounds
    loint = min(intervals | set([modint]), key=lambda y: y[1])
    hiint = max(intervals | set([modint]))
    l_idx = get_indices(xs, loint)
    h_idx = get_indices(xs, hiint)

    # recurse low
    if dip.diptst(xs[:l_idx[1]+1], NTRIALS)[1] < alpha:
        intervals |= unidip(xs[:l_idx[0]+1], alpha, False, plotdat=plotdat)

    # recurse high
    if dip.diptst(xs[h_idx[0]:], NTRIALS)[1] < alpha:
        intervals |= unidip(xs[h_idx[1]:], alpha, False, plotdat=plotdat)

    # return all intervals
    return intervals


def get_indices(arr, vals):
    """ Get indices from numpy array (arr) from values (vals) """
    sorter = np.argsort(arr)
    return sorter[np.searchsorted(arr, vals, sorter=sorter)]


def _db_plot(dat, sub, ints):
    """ Plot complete data, highlight subset currently being searched,
    and add vertical lines for discovered intervals. (only intervals of
    the current level appear.)
    """
    import matplotlib.pyplot as plt
    plt.style.use("seaborn")
    axs = plt.hist(dat, bins=30, density=True)
    plt.axvspan(sub[0], sub[1], color="orange", alpha=.3)
    plt.vlines(list(ints), 0, max(axs[0]))
    plt.show()


def _get_full_interval(dat, modint):
    """ Expands discovered intervals
    When looking at unimodal data the dip test tends to return a very narrow
    interval, which can lead to conflicts later. This tends to happen after
    recursing left or right.

    Our solution, taken from the original unidip, is to mirror the data such
    that it is bimodal. We are then able to fully capture the mode, and return
    the full mode from the original data.
    """
    ldat = _mirror_data(dat, left=True)
    ldip = dip.diptst(ldat, NTRIALS)
    rdat = _mirror_data(dat)
    rdip = dip.diptst(rdat, NTRIALS)

    if ldip[0] > rdip[0]:
        full_indxs = _un_mirror_idxs(get_indices(ldat, ldip[2]),
                                     len(dat), modint, True)
    else:
        full_indxs = _un_mirror_idxs(get_indices(rdat, rdip[2]),
                                     len(dat), modint, False)
    return tuple(dat[full_indxs])


def _mirror_data(dat, left=False):
    """ Mirror dataset
    input: [1, 2, 3] output: [-2, -1, 0, 1, 2]
    """
    wdat = np.array(dat)
    if left:
        pivot = np.min(wdat)
        sdat = wdat-pivot
        mdat = np.concatenate((-sdat[sdat > 0], sdat))
    else:
        pivot = np.max(wdat)
        sdat = wdat-pivot
        mdat = np.concatenate((sdat, -sdat[sdat < 0]))
    return np.sort(mdat)


def _un_mirror_idxs(vals, length, modint, left):
    """ wrapper for _un_mirror_idx() """
    if vals[0] < length and vals[1] > length:
        return modint

    ori_idxs = np.array([_un_mirror_idx(i, length, left) for i in vals])
    return ori_idxs if ori_idxs[0] < ori_idxs[1] else np.flip(ori_idxs, -1)


def _un_mirror_idx(idx, length, left):
    """ take index from mirrored data and return the appropriate index for
    original data.

    Using indices prevents issues with floating point imprecision.
    """
    if left:
        return idx-length if idx >= length else length-idx
    else:
        return (2 * length)-idx if idx > length else idx


def _test(filename, plot=False, debug=False, **kwargs):
    """ test filename's.csv peakitude :) """
    from time import time
    dat = np.genfromtxt(filename, delimiter=",")
    print(f"length of test on {filename}: {len(dat)}")
    start = time()
    if debug:
        ints = unidip(dat, plotdat=dat, **kwargs)
    else:
        ints = unidip(dat, **kwargs)
    end = time()
    print(f"# intervals returned: {len(ints)}")
    for i in sorted(ints):
        print(i)
    print(f"time taken {end-start:.2f}sec\n")
    if plot:
        import matplotlib.pyplot as plt
        ax = plt.hist(dat, bins=30)
        plt.vlines(list(ints), 0, max(ax[0]))
        plt.show()


if __name__ == "__main__":
    _test("tests/testsmall.csv")
    _test("tests/peak1.csv", False)
    _test("tests/peak2.csv", False)
    _test("tests/peak3.csv", False)
    _test("tests/large3.csv", False)
    _test("tests/test10p.csv")
    _test("tests/test1or10p.csv", alpha=.3)
    _test("tests/test0.5sig.csv")
    print("finished testing!")
