from GammaAP import *
#apfilename
#ranal

if __name__ == '__main__':
	ga = GammaAP()
	ga.fullAnalysisLoadAP4(sys.argv[1], analyzevm=-1, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800, vmnoise=1)
	
