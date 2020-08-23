from GammaAP import *
ap = GammaAP()
ap.loadnormalizedAP4("RES1_TBS3600_R2_DQ1_EB100-10000.ap.sim.ap.ap4.compare")
ap.calculateLSexp(verbose=1, plot=1, rescol=3)
