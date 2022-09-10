import os,sys,math,random
import numpy as np

start = int(sys.argv[1])
nruns = int(sys.argv[2])

for i in range(start, start+nruns, 1):
	
	dirname = str(i).zfill(10)
	os.system("mkdir " + dirname)
	os.chdir(dirname)
	os.system("cp ../BINS/BIN? .")
	os.system("cp ../sources.multi .")
	os.system("python ../generate_lc.py 1 .")
	for i in range(0, 7):
		j = "BIN" + str(i + 1)
		os.system("map.rb FM3.119_ASDCSTD_H0025 "+j+"R 0 0 355.447 -0.2689 emin=100 emax=10000 dq=1 eb=0 binsize=0.1 timelist="+j+" useEDPmatrixforEXP=0 skytype=4")
		os.system("multi.rb FM3.119_ASDCSTD_H0025 "+j+"R.maplist4 sources.multi "+j+"RES galmode2=3 isomode2=3 ulcl=1 fluxcorrection=1")
	os.chdir("..")

