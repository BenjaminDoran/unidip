import unidip.unidip as unidip


unidip._test("./tests/testsmall.csv") # 0 peaks, not enough data
unidip._test("./tests/peak1.csv") # 1 peak
unidip._test("./tests/peak2.csv", debug=False) # 2 peaks
unidip._test("./tests/peak3.csv", debug=False) # 3 peaks
unidip._test("./tests/large3.csv", debug=False) # 3 peaks
unidip._test("./tests/test10p.csv", plot=False) # 10 peaks
unidip._test("./tests/test1or10p.csv", plot=False, alpha=0.3) # 10 peaks small gaps
unidip._test("./tests/test0.5sig.csv", debug=False, alpha=.05) # 5 peaks
unidip._test("./tests/histnotsig.csv", is_hist=True) # 3 peaks, but n of 10, so 1
unidip._test("./tests/hist3p.csv", debug=False, is_hist=True) # 3 peaks, n of 60
# unidip._test("./tests/negEntIdxErr.csv", is_hist=True, debug=True) # off by one error in diptst / fixed
# unidip._test("./tests/negEntMaxRecErr.csv", is_hist=True, debug=True) # lead to errs because search int finds ~3 positions of mode on each side of inter-modal space
print("finished testing!")