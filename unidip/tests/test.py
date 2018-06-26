from unittest import TestCase

import unidip

EXAMPLES = "./unidip/tests/example_data/"

class Tester(TestCase):
    # numeric tests
    def test_0_peaks(self):
        """ testing 0 peaks, not enough data """
        pn = 0
        ints = unidip.test_unidip(f"{EXAMPLES}testsmall.csv") 
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    def test_1_peak(self):
        """ test 1 peak"""
        pn = 1
        ints = unidip.test_unidip(f"{EXAMPLES}peak1.csv", plot=False)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    def test_2_peaks(self):
        """ test 2 peaks """
        pn = 2
        ints = unidip.test_unidip(f"{EXAMPLES}peak2.csv", plot=False, debug=False)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    def test_3_peaks(self): 
        """ test 3 peaks """
        pn = 3
        ints = unidip.test_unidip(f"{EXAMPLES}peak3.csv", plot=False)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    def test_large3_peaks(self):
        """ test large three peaks """
        pn = 3
        ints = unidip.test_unidip(f"{EXAMPLES}large3.csv", plot=False)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    def test_10_peaks(self):
        """ test 10 peaks """
        pn = 10
        ints = unidip.test_unidip(f"{EXAMPLES}test10p.csv", plot=False, debug=False)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    def test_1or10_peaks(self):
        """ test hidden 10 peaks """
        pn = 10
        ints = unidip.test_unidip(f"{EXAMPLES}test1or10p.csv", plot=False, alpha=0.32, mrg_dst=5)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    def test_5wide_peaks(self):
        """ test 5 wide peaks """
        pn = 5
        ints = unidip.test_unidip(f"{EXAMPLES}test0.5sig.csv", plot=False, debug=False, alpha=.05)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    # histogram tests
    def test_histsmall_peak(self):
        """ Tests on 10 density points from histogram """
        pn = 1
        ints = unidip.test_unidip(f"{EXAMPLES}histnotsig.csv", plot=False, is_hist=True)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
    def test_hist3_peaks(self):
        """ test 3 histogram peaks """
        pn = 3
        ints = unidip.test_unidip(f"{EXAMPLES}hist3p.csv", plot=False, is_hist=True)
        assert len(ints) == pn, f"Found wrong number of peaks! Should be {pn}, is {len(ints)}."
