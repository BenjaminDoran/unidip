"""
Author: Benjamin Doran
Date: Jan 2018

Algorithm ref:
"Skinny-dip: Clustering in a sea of noise" by Samuel Maurus and Claudia Plant.
http://www.kdd.org/kdd2016/subtopic/view/skinny-dip-clustering-in-a-sea-of-noise
"""

import numpy as np
import dip

NTRIALS = 100
MERGE_DISTANCE = 1

# wrapper function in case further error checking is desired
def unidip(dat, is_hist=False, alpha=0.05, numt=NTRIALS, plotdat=None):
    """ Perform unidip algorithm on 1d array

        INPUT:
        idxs: x axis of data, default for data like that created by
                np.random.normal() or np.random.uniform()
        hist: density along x axis, allows for dipping across discrete functions
                such as a poisson distribution.
        alpha: tuning parameter, sets significance level of p_values
        plotdat: determines whether to plot the data at each recursion level,
                    set to same as dat argument

        RETURNS:
        set of tuples: each tuple containing the start and end positions on the
                        x axis.
    """
    offset = 0
    is_model = True
    _dat_is_x = not is_hist
    modidxs = _unidip(dat, offset,
                      _dat_is_x, alpha, is_model, numt, plotdat)
    return merge_intervals(modidxs)

def merge_intervals(idxs, merge_distance=MERGE_DISTANCE):
    """ merge intervals that are touching """
    midxs = []
    for idx in sorted(idxs):
        if not midxs:
            midxs.append(idx)
        else:
            lower = midxs[-1]
            # adjacent or overlapping (adjust MERGE_DISTANCE to merge intervals)
            # that are close enough
            if idx[0] - lower[1] <= merge_distance:
                midxs[-1] = (lower[0], idx[1])
            else:
                midxs.append(idx)
    return midxs


def _unidip(dat, offset=0, _dat_is_x=True, alpha=0.05,
            _is_model=True, numt=NTRIALS, plotdat=None):
    """ Perform unidip algorithm on 1d array

        INPUT:
        dat: x axis of data, default for data created by
                np.random.normal() or np.random.uniform() and the like
        idxs: density along x axis, allows for dipping across discrete functions
                such as a poisson distribution by inputting the tally of bin
                heights
        alpha: tuning parameter, sets significance level of p_values
        _is_model: internal should not be changed
        plotdat: determines whether to plot the data at each recursion level,
                set to same as dat argument

        RETURNS:
        set of tuples: each tuple containing the start and end positions on the
                        x axis.
    """
    # _idxs = np.arange(len(dat))
    if _dat_is_x:
        dat = np.msort(dat)
    interval_idxs = list()
    _, pval, modidx, modint = dip.diptst(dat, _dat_is_x, numt)

    if not plotdat is None: # if plotting -> show intervals
        if _dat_is_x:
            _db_plot(plotdat, (dat[0], dat[-1]), [modint], _dat_is_x)
        else:
            _db_plot(plotdat, (offset, offset+len(dat)-1),
                     [(offset+modidx[0], offset+modidx[1])], _dat_is_x)

    # not enough data to count it as significant
    if pval is None:
        return []
    # is unimodal, return interval
    elif pval > alpha:
        if _is_model:
            interval_idxs.append((offset, offset+len(dat)-1))
        else:
            wideint, wideidx = _get_full_interval(dat, modint, _dat_is_x, numt)
            interval_idxs.append((offset+wideidx[0], offset+wideidx[1]))
        return interval_idxs

    # recurse into model interval
    rmidx = _unidip(dat[modidx[0]:modidx[1]+1], offset+modidx[0],
                    _dat_is_x, alpha, True, numt, plotdat)

    # add returned intervals to our collection
    interval_idxs += rmidx
    # undo offset to get correct indices to data in recursion layer
    subd = list(map(lambda t: (t[0]-offset, t[1]-offset), interval_idxs))
    # upper and lower bounds
    l_idx = min(subd + [modidx], key=lambda y: y[1])
    h_idx = max(subd + [modidx])

    # recurse low
    pval = dip.diptst(dat[:l_idx[1]+1], _dat_is_x, numt)[1]
    if not pval is None and pval < alpha:
        rlidx = _unidip(dat[:l_idx[0]+1], offset,
                        _dat_is_x, alpha, False, numt, plotdat)
        interval_idxs += rlidx

    # recurse high
    pval = dip.diptst(dat[h_idx[0]:], _dat_is_x, numt)[1]
    if not pval is None and pval < alpha:
        rhidx = _unidip(dat[h_idx[1]:], offset+h_idx[1],
                        _dat_is_x, alpha, False, numt, plotdat)
        interval_idxs += rhidx

    # return all intervals
    return interval_idxs


def _db_plot(dat, sub, ints, is_idxs=True):
    """ Plot complete data, highlight subset currently being searched,
        and add vertical lines for discovered intervals. (only intervals of
        the current level appear.)
    """
    import matplotlib.pyplot as plt
    plt.style.use('seaborn')
    if is_idxs:
        plt.hist(dat, bins=30)
    else:
        plt.step(list(range(len(dat))), dat)
        plt.fill_between(list(range(len(dat))), dat, step="pre", alpha=.4)
    plt.axvspan(sub[0], sub[1], color="orange", alpha=.3)
    for i in ints:
        plt.axvspan(i[0], i[1], color="green", alpha=.1)
    for i in ints:
        plt.axvline(i[0], color="black")
        plt.axvline(i[1], color="black")
    plt.show()


def _get_full_interval(dat, modint, is_idxs, numt):
    """ Expands discovered intervals
        When looking at unimodal data the dip test tends to return a very narrow
        interval, which can lead to conflicts later. This tends to happen after
        recursing left or right.

        Our solution, taken from the original unidip, is to mirror the data such
        that it is bimodal. We are then able to fully capture the mode, and return
        the full mode from the original data.
    """
    ldat = _mirror_data(dat, is_idxs, left=True)
    ldip = dip.diptst(ldat, is_idxs, numt)
    rdat = _mirror_data(dat, is_idxs, left=False)
    rdip = dip.diptst(rdat, is_idxs, numt)

    if ldip[0] > rdip[0]:
        full_indxs = _un_mirror_idxs(ldip[2], len(dat), modint, True)
    else:
        full_indxs = _un_mirror_idxs(rdip[2], len(dat), modint, False)
    return (tuple(dat[full_indxs]), tuple(full_indxs))


def _mirror_data(dat, is_idxs=True, left=False):
    """ Mirror dataset
    input: [1, 2, 3] output: [-2, -1, 0, 1, 2]
    """
    wdat = np.array(dat)
    if is_idxs:
        if left:
            pivot = np.min(wdat)
            sdat = wdat-pivot
            mdat = np.concatenate((-sdat[sdat > 0], sdat))
            np.sort(mdat)
        else:
            pivot = np.max(wdat)
            sdat = wdat-pivot
            mdat = np.concatenate((sdat, -sdat[sdat < 0]))
            np.sort(mdat)
    else:
        if left:
            mdat = np.concatenate((np.flip(wdat[1:], -1), wdat))
        else:
            mdat = np.concatenate((wdat[:-1], np.flip(wdat, -1)))
    return mdat


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
        idxs = idx-length if idx >= length else length-idx
    else:
        idxs = (2 * length)-idx if idx > length else idx
    return idxs


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
    print("interval idxs:")
    for i in sorted(ints):
        print(i)
    print(f"time taken {end-start:.2f}sec\n")
    if plot:
        import matplotlib.pyplot as plt
        plt.style.use('seaborn')
        plt.hist(dat, bins=30)
        for i in ints[0]:
            plt.axvspan(i[0], i[1], color="green", alpha=.3)
        for i in ints[0]:
            plt.axvline(i[0])
            plt.axvline(i[1])
        plt.show()


if __name__ == "__main__":
    # _test("tests/testsmall.csv") # 0 peaks, not enough data
    # _test("tests/peak1.csv") # 1 peak
    # _test("tests/peak2.csv", debug=False) # 2 peaks
    # _test("tests/peak3.csv", debug=False) # 3 peaks
    # _test("tests/large3.csv", debug=False) # 3 peaks
    # _test("tests/test10p.csv", plot=False) # 10 peaks
    _test("tests/test1or10p.csv", plot=False, alpha=0.3) # 10 peaks small gaps
    # _test("tests/test0.5sig.csv", debug=False, alpha=.05) # 5 peaks
    # _test("tests/histnotsig.csv", is_hist=True) # 3 peaks, but n of 10, so 1
    # _test("tests/hist3p.csv", debug=False, is_hist=True) # 3 peaks, n of 60
    # _test("tests/negEntIdxErr.csv", is_hist=True, debug=True)
    print("finished testing!")
