from GammaAP import *

#1 ap file
#2 ranal

if __name__ == '__main__':
	ap = GammaAP()
	print(sys.argv[1] + " " + str(sys.argv[2]) + " " + str(0.00054))
	ap.normalizeAP3(sys.argv[1], sys.argv[2], 0.00054)
