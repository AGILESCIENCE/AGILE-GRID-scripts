from GammaAP import *

if __name__ == '__main__':
	ga = GammaAP()
	ga.fullAnalysis3(sys.argv[1], 0, 48, 0.5e-06, 5.0e-06, 100, 10000, 10800)
