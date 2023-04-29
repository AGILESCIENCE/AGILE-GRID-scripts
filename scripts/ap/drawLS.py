from GammaAP import *

def main():
	ap = GammaAP()
	ap.evaluateLS(sys.argv[1], rescol=int(sys.argv[2]), plot=1)

if __name__ == '__main__':
	main()