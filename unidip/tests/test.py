from unittest import TestCase

import unidip

EXAMPLES = "./unidip/tests/example_data/"

class Tester(TestCase):
    # numeric tests
    def test_0_peaks(self):
        """ testing 0 peaks, not enough data """
        unidip.test_unidip(f"{EXAMPLES}testsmall.csv") 
    def test_1_peak(self):
        """ test 1 peak"""
        unidip.test_unidip(f"{EXAMPLES}peak1.csv", plot=False)
    def test_2_peaks(self):
        """ test 2 peaks """
        unidip.test_unidip(f"{EXAMPLES}peak2.csv", plot=False, debug=False)
    def test_3_peaks(self): 
        """ test 3 peaks """
        unidip.test_unidip(f"{EXAMPLES}peak3.csv", plot=False)
    def test_large3_peaks(self):
        """ test large three peaks """
        unidip.test_unidip(f"{EXAMPLES}large3.csv", plot=False)
    def test_10_peaks(self):
        """ test 10 peaks """
        unidip.test_unidip(f"{EXAMPLES}test10p.csv", plot=False, debug=False)
    def test_1or10_peaks(self):
        """ test hidden 10 peaks """
        unidip.test_unidip(f"{EXAMPLES}test1or10p.csv", 
                           plot=False, alpha=0.3, merge_distance=5)
    def test_5wide_peaks(self):
        """ test 5 wide peaks """
        unidip.test_unidip(f"{EXAMPLES}test0.5sig.csv",
                           plot=False, debug=False, alpha=.05)
    # histogram tests
    def test_histsmall_peak(self):
        """ Tests on 10 density points from histogram """
        unidip.test_unidip(f"{EXAMPLES}histnotsig.csv", 
                          plot=False, is_hist=True)
    def test_hist3_peaks(self):
        """ test 3 histogram peaks """
        unidip.test_unidip(f"{EXAMPLES}hist3p.csv", plot=False, is_hist=True)
    # Neg Entropy data
    def test_negEntIdxErr_peaks(self):
        """ tests finding 5 hist peaks """
        unidip.test_unidip(f"{EXAMPLES}negEntIdxErr.csv", 
                           is_hist=True, plot=False)
    def test_negEntMaxRecErr(self):
        """ Test for max recursion error on 4 peaks"""
        unidip.test_unidip(f"{EXAMPLES}negEntMaxRecErr.csv", 
                           is_hist=True, plot=False)