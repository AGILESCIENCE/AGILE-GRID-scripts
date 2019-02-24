from GammaAP import *

if __name__ == '__main__':
	ap = GammaAP()
	ap.plotVonMisses(sys.argv[1], float(sys.argv[2]), 1, 2)
