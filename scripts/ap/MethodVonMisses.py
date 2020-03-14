# DESCRIPTION
#       Agileap: AGILE Observatory Aperture Photometry Analysis
# NOTICE
#      Any information contained in this software
#      is property of the AGILE TEAM and is strictly
#      private and confidential.
#      Copyright (C) 2005-2020 AGILE Team.
#          Bulgarelli Andrea <andrea.bulgarelli@inaf.it>
#          Valentina Fioretti <valentina.fioretti@inaf.it>
#          Parmiggiani Nicolò <nicolo.parmiggiani@inaf.it>
#      All rights reserved.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import string, os, sys
import numpy as np

class MethodVonMisses:
	
	def __init__(self, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100):
		self.apfile = ""
		self.freqmin = float(freqmin)
		self.freqmax = float(freqmax)
		self.vmnumax = float(vmnumax)
		return
	
	def plotVonMisses(self, filename, mu=0.0, verbose=1, plot=1):
		diml=0
		with open(filename, "r") as ins:
			for line in ins:
				if(line != ""):
					val = line.split()
					#print(len(val))
					if len(val) > 2:
						muv = float(val[0])
						if(muv == float(mu)):
							diml += 1
		
		#print("diml " + str(diml))
		frequency = np.zeros(diml)
		power = np.zeros(diml)
		
		diml=0
		
		maxf = 0
		maxp = 0
		with open(filename, "r") as ins:
			for line in ins:
				if(line != ""):
					val = line.split()
					if len(val) > 2:
						#print(str(val[0]) + " " + str(val[1]) + " " + str(val[2]))
						muv = float(val[0])
						if(muv == mu):
							frequency[diml] = float(val[1])
							power[diml] = float(val[2])
							
							if maxf == 0:
								maxf = frequency[diml]
								maxp = power[diml]
							
							if power[diml] > maxp:
								maxf = frequency[diml]
								maxp = power[diml]
							
							diml += 1
		
		daymax = 1 / maxf / 86400.
		
		if verbose == 1:
			print('pls -1 with power of ' + str(maxp) + ' at frequency ' + str(maxf) + ' Hz (' + str(daymax) + ') days' )
		
		if plot > 0:
			plt.plot(frequency, power)
			plt.gca().set_xscale("log")
			#plt.gca().set_yscale("log")
			if plot == 1:
				plt.show()
			if plot == 2:
				plt.savefig(filename+"."+str(mu)+".png", format='png', dpi=1200)
				
		return maxf, maxp

	def scanVM(self, vmfilename, coltype):
		vmtype = "vm" + str(coltype)
		fileclean = open(self.apfile + ".apres", "a")
		muA = []
		mu = 0.0
		muA.append(mu)
		with open(vmfilename, "r") as ins:
			for line in ins:
				if(line != ""):
					val = line.split()
					#print(len(val))
					if len(val) > 2:
						muv = float(val[0])
						if(muv != mu):
							muA.append(muv)
							mu = muv
		
		for i in muA:
			maxf, maxp = self.plotVonMisses(vmfilename, i, 0, 0)
			daymax = 1 / maxf / 86400.
			print(vmtype + 'mu ' + str(i) + ' -1 ' + str(maxp) + ' ' + str(maxf) + ' ' + str(daymax))
			fileclean.write(vmtype + 'mu ' + str(i) + ' -1 ' + str(maxp) + ' ' + str(maxf) + ' ' + str(daymax)+ "\n")
			fileclean.close()
		return
	
	
	def runVomMisses(self, nthreads, ii):
		self.vmnoise=1
		print(self.apfile)
		#Arguments: threads number, noise scalability flag (zero for NOT to scale), f_min, f_max, nu_min, nu_max
		cmd = "module load icc-18.0.1; module load gcc-5.4.0; "+os.environ['AGILE']+"/bin/eval_vonmises.prg "+str(nthreads)+" "+ str(self.vmnoise) + " " + str(self.freqmin) + " " + str(self.freqmax) + " 0 " + str(self.vmnumax) + " < " + self.apfile + ".vm"+str(ii)+" > " + self.apfile + ".vm"+str(ii)+".res"
		print(cmd)
		os.system(cmd)
	
	
	def runVomMissesGridFreq(self, ii, ngrid=1000, tspan=10800):
		#Arguments: f_min, f_max, N, original time series span T
		#T is used to organize an almost-logarithmic f-scale; use any T<=0 to request linear f-scale
		cmd = "module load icc-18.0.1; module load gcc-5.4.0; "+os.environ['AGILE']+"/bin/grid_freq.prg "+ str(self.freqmin) + " " + str(self.freqmax) + " " + str(ngrid) + " " + str(tspan) + " < " + self.apfile + ".vm"+str(ii)+".res > " + self.apfile + ".vm"+str(ii)+".resgf"
		print(cmd)
		os.system(cmd)
	
	
	def significanceVonMisses(self, nthreads, ii, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100):
		print("significance von misses periodogram")
		cmd = "module load icc-18.0.1; module load gcc-5.4.0; "+os.environ['AGILE']+"/bin/coeffs_XY.prg "+str(nthreads)+ " " + str(freqmin) + " " + str(freqmax) + " 0 " + str(vmnumax) + " < " + self.apfile + ".vm"+str(ii)+" > " + self.apfile + ".vm"+str(ii)+".sig"
		print(cmd)
		os.system(cmd)
	
	
	def evaluateSignificanceVonMisses(self, W, Xnumax, Ynumax, z, Ynu0=1):
		#z max  value of the periodogram
		
		#res =W * pow(2.7182818, -z) * (z * 1.20555 + (5.116071 + 1) * (sqrt(z) / 2) )
		#W * e^(-z) * (  z * Xnumax + ( Ynumax + Ynu0 ) * (sqrt(z) / 2)  )
		#import numpy as np
		#W =
		#Xnumax =
		#Ynumax =
		#Ynu0=
		#z=
		#W * np.power(np.e, -z) * (  z * Xnumax + ( Ynumax + Ynu0) * np.sqrt(z) / 2.0 )
		
		#import numpy as np
		#W =
		#Xnumax =
		#Ynumax =
		#z=
		#W * np.power(np.e, -z) * (  2 * z * Xnumax + Ynumax  * np.sqrt(z) )
		#W * e^(-z) * (  2 * z * Xnumax + Ynumax  * sqrt(z) )
		#res2=W * pow(2.7182818, -z) * (2 * z * 1.20555 + 5.116071 * sqrt(z))
		sig1=W * np.power(np.e, -z) * (  z * Xnumax + ( Ynumax + Ynu0) * np.sqrt(z) / 2.0 )
		print("sig1= " + str(sig1))
		sig2=W * np.power(np.e, -z) * (  2 * z * Xnumax + Ynumax  * np.sqrt(z) )
		print("sig2= " + str(sig2))

	def runVomMissesStep(self, ii, vonmissesthread=48, ngridfreq=1000, tgridfreq=10800):
		self.runVomMisses(vonmissesthread, ii)
		self.runVomMissesGridFreq(ii, ngridfreq, tgridfreq)
		self.scanVM(self.apfile + ".vm"+str(ii)+".resgf", ii)
	
	def fullAnalysis(self, apfile, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800):
		
		self.apfile = apfile
		self.freqmin = float(freqmin)
		self.freqmax = float(freqmax)
		self.vmnumax = float(vmnumax)
		
		self.runVomMissesStep(0, vonmissesthread, ngridfreq, tgridfreq)
		self.runVomMissesStep(1, vonmissesthread, ngridfreq, tgridfreq)
		self.runVomMissesStep(2, vonmissesthread, ngridfreq, tgridfreq)
		self.runVomMissesStep(3, vonmissesthread, ngridfreq, tgridfreq)
		self.runVomMissesStep(8, vonmissesthread, ngridfreq, tgridfreq)
		

