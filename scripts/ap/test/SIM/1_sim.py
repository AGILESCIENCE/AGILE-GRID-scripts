import os
from SimAP import *
sa = SimAP()
sa.calculateAPLC(mode=0, apfile="./RES1_TBS3600_R2_DQ1_EB100-10000.ap", ranal=2, gasvalue=0.00054, gal=0.5, iso=6, emin=100, emax=10000, period = 729855.36, peak_size=300e-08, gindex=2.2)
sa.calculateAPLC(mode=0, apfile="./RES1_TBS14400_R2_DQ1_EB100-10000.ap", ranal=2, gasvalue=0.00054, gal=0.5, iso=6, emin=100, emax=10000, period = 729855.36, peak_size=300e-08, gindex=2.2)

#print('compare RES1_TBS14400_R2_DQ1_EB100-10000.ap.sim.ap')
#os.system("diff RES1_TBS14400_R2_DQ1_EB100-10000.ap.sim.ap RES1_TBS14400_R2_DQ1_EB100-10000.ap.sim.ap.compare")
#print('compare RES1_TBS3600_R2_DQ1_EB100-10000.ap.sim.ap')
#os.system("diff RES1_TBS3600_R2_DQ1_EB100-10000.ap.sim.ap RES1_TBS3600_R2_DQ1_EB100-10000.ap.sim.ap.compare")
