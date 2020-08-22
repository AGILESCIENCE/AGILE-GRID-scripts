import os

from GammaAP import *
ap = GammaAP()
ap.normalizeAP("E1dq1_86400s_emin100_emax10000_r4.ap", ranal=4, gasvalue=0.0, gal=0.0, iso=10, gindex=2.1, evalULalgorithm=1)

print("check compare")
os.system("diff E1dq1_86400s_emin100_emax10000_r4.ap.ap4 E1dq1_86400s_emin100_emax10000_r4.ap.ap4.compare")
