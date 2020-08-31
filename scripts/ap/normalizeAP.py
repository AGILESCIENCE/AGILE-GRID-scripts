from GammaAP import *

#1 ap file
#2 ranal

if __name__ == '__main__':
	ap = GammaAP()

	ap.normalizeAP(sys.argv[1], ranal=sys.argv[2], gasvalue=0.00054, gal=0.7, iso=10, emin=100, emax=10000, gindex=2.1, writevonmissesfiles=0, evalULalgorithm=1)
