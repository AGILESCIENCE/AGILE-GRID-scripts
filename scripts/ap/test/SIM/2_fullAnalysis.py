import os
from GammaAP import *
ap = GammaAP()
#rateMeanBkgExpected=0.000001301365885870
ap.fullAnalysis("RES1_TBS14400_R2_DQ1_EB100-10000.ap.sim.ap.compare", analyzevm=1, rateMeanBkgExpected=-1, ranal=2, gasvalue=0.00056, gal=0.5, iso=6, gindex=2.2, writevonmissesfiles=1, tgridfreq=14400, freqmin=1.0e-07, freqmax=5.0e-06, vmnumax=100, ngridfreq=10000, vmnoise=1)

print("diff...")

os.system("diff RES1_TBS14400_R2_DQ1_EB100-10000.ap.sim.ap.compare.ap4  RES1_TBS14400_R2_DQ1_EB100-10000.ap.sim.ap.compare.ap4.compare")
