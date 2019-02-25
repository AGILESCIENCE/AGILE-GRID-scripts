from GammaAP import *

if __name__ == '__main__':
	ap = GammaAP()
	#last parameter
	#0 do not draw
	#1 draw on screen
	#2 draw on file
	ap.plotVonMisses(sys.argv[1], float(sys.argv[2]), 1, 2)
