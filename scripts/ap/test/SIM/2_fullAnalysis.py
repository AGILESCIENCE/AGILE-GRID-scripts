from GammaAP import *
ap = GammaAP()
ap.fullAnalysis("RES1_TBS14400_R2_DQ1_EB100-10000.ap.sim.ap", ranal=2, gasvalue=0.00056, gal=0.5, iso=6, gindex=2.2, writevonmissesfiles=1, tgridfreq=14400, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000)
