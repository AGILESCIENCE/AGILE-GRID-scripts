from GammaAP import *

if __name__ == '__main__':
	ap = GammaAP()
	ap.evaluateLS(sys.argv[1], rescol=int(sys.argv[2]), plot=1)
