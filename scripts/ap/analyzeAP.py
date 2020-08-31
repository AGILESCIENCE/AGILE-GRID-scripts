from GammaAP import *
#apfilename
#ranal

if __name__ == '__main__':
	ga = GammaAP()
	ga.fullAnalysis(sys.argv[1], ranal=sys.argv[2], gasvalue=0.00054, analyzevm=0, writevonmissesfiles=0, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=10000, tgridfreq=10800, vmnoise=1, gal=0.7, iso=10, emin=100, emax=10000, gindex=2.1,  evalULalgorithm=1)
):