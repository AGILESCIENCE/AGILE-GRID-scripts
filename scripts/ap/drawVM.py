from GammaAP import *

def main():
	ap = GammaAP()
	#last parameter
	#0 do not draw
	#1 draw on screen
	#2 draw on file
	ap.plotVonMisses(sys.argv[1], float(sys.argv[2]), 0, 2)

if __name__ == '__main__':
	main()