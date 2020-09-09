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
#          Alessio Aboudan
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
import math
import matplotlib.pyplot as plt
from EvalRates import *
from MethodVonMisses import *
from MethodLS import *
from APSignificance import *
from MathUtils import *

#gasvalue IGR = 0.00054
#gasvalue Vela = 0.00018
#TODO
#1) fare un calcolo del ratewmean basato su una media mobile

class NormalizeAP:

	def __init__(self):
		return

	# Section 2
	#
	# [1] Robin Corbet, LAT Light Curve Analysis Aperture Photometry and Periodicity Searches, slides at COSPAR CBW 2010
	# [2] Gehrels, 1986, ApJ, 303, 336
	# [3] Kraft, Burrows, & Nousek, 1991, ApJ, 374, 344
	# Variance weighting corrected for low count statistic - Weighted Power Spectrum	 [1] [2]

	def normalizeAB3(self, expdataA, ctsdataA, rateBkgExpected, normrateAB11, normrateAB12, normrateAB13, normrateAB14, normrateAB21, normrateAB22, normrateAB23, normrateAB24, normrateAB11aa, normrateAB21aa, ratediffR1, ratediffR2, ratediffR3, ratediffR4, ratediffR1AA, rate, rateError, expBasedRateError):
		sum1 = 0.0
		sum2 = 0.0
		sum3 = 0.0
		sum4 = 0.0
		sum5 = 0.0
		diml = len(expdataA)
		#rate = np.zeros(diml)
		#rateError = np.zeros(diml)
		
		n=0
		for e_i in expdataA:
			
			n_i = ctsdataA[n] #cts
			
			rate_i = n_i / e_i #cts / cm2 s
			rate[n] = rate_i #cts / cm2 s
			
			s_ip = 0.5 + np.sqrt(n_i + 0.25)
			s_im = -0.5 + np.sqrt(n_i + 0.25)
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i
			rateError[n] = s_i

			sum1 += rate_i / (s_i*s_i)
			sum2 += 1.0 / (s_i*s_i)
			sum3 += rate_i
			sum4 += n_i
			sum5 += e_i
			
			n = n + 1
			
		rateWeightedMean1 = sum1 / sum2 #cts / cm2 s
		rateWeightedMean2 = sum3 / n #cts / cm2 s
		rateWeightedMean3 = sum4 / sum5 #cts / cm2 s
		rateWeightedMean4 = rateBkgExpected #cts / cm2 s
		print("rateWeightedMean1:          %.3e"% rateWeightedMean1)
		print("rateWeightedMean2:          %.3e"% rateWeightedMean2)
		print("rateWeightedMean3:          %.3e"% rateWeightedMean3)
		print("rateWeightedMean4 bkgmodel: %.3e"% rateWeightedMean4)
		print("rateWeightedMean4 bkgmodel: %.19f"% rateWeightedMean4)
		
		n=0
		sum1 = 0.0
		sum2 = 0.0
		for rate_i in rate:
			e_i = expdataA[n]
			ctspred_i = e_i * rateWeightedMean1 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			sum1 += rate_i / (sp_i*sp_i)
			sum2 += 1.0 / (sp_i*sp_i)
		
		rateWeightedMean1aa = sum1 / sum2 #cts / cm2 s
		print("rateWeightedMean1aa:        %.3e" %rateWeightedMean1aa)
		

		n=0
		for rate_i in rate:
			e_i = expdataA[n]
			s_i = rateError[n]
			
			###############7.1.1
			#res column 0 normrateAB11
			ctspred_i = e_i * rateWeightedMean1 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB11[n] = (rate_i - rateWeightedMean1) / (sp_i*sp_i)
			
			###############7.1.2
			#res column 1 normrateAB12
			ctspred_i = e_i * rateWeightedMean2 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB12[n] = (rate_i - rateWeightedMean2) / (sp_i*sp_i)
			
			###############7.1.3
			#res column 2 normrateAB13
			ctspred_i = e_i * rateWeightedMean3 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			expBasedRateError[n] = sp_i
			normrateAB13[n] = (rate_i - rateWeightedMean3) / (sp_i*sp_i)
			
			###############7.1.4
			#res column 3 normrateAB14
			ctspred_i = e_i * rateWeightedMean4 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB14[n] = (rate_i - rateWeightedMean4) / (sp_i*sp_i)
			
			###############7.2.1
			#res column 4 normrateAB21
			ctspred_i = e_i * rateWeightedMean1 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB21[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.2
			#res column 5 normrateAB22
			ctspred_i = e_i * rateWeightedMean2 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB22[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.3
			#res column 6 normrateAB23
			ctspred_i = e_i * rateWeightedMean3 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB23[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.4
			#res column 7 normrateAB24
			ctspred_i = e_i * rateWeightedMean4 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB24[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.3.1
			#rateAB31[n] = (rate_i - rateWeightedMean1) / (s_i*s_i)
			
			###############7.3.2
			#rateAB32[n] = (rate_i - rateWeightedMean2) / (s_i*s_i)
			
			###############7.3.3
			#rateAB33[n] = (rate_i - rateWeightedMean3) / (s_i*s_i)
			
			###############7.3.4
			#rateAB34[n] = (rate_i - rateWeightedMean4) / (s_i*s_i)
			
			###############7.1.1aa
			#res column 8 normrateAB11aa
			ctspred_i = e_i * rateWeightedMean1aa #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB11aa[n] = (rate_i - rateWeightedMean1aa) / (sp_i*sp_i)
			
			###############7.2.1aa
			#res column 9 normrateAB21aa
			ctspred_i = e_i * rateWeightedMean1aa #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			normrateAB21aa[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.4.1
			#rateAB41[n] = (rate_i ) / (s_i*s_i)
			
			###############rate
			#res column 10
			ratediffR1[n] = (rate_i - rateWeightedMean1) 
			
			###############rate
			#res column 11
			ratediffR2[n] = (rate_i - rateWeightedMean2)
			
			###############rate
			#res column 12
			ratediffR3[n] = (rate_i - rateWeightedMean3)
			
			###############rate
			#res column 13
			ratediffR4[n] = (rate_i - rateWeightedMean4)
			
			#res column 14
			ratediffR1AA[n] = (rate_i - rateWeightedMean1aa)
			
			n = n + 1
		
		return rateWeightedMean1, rateWeightedMean2, rateWeightedMean3, rateWeightedMean4, rateWeightedMean1aa


#########################################
class GammaAP:

	def __init__(self):
		self.ncols = 48
		self.diml = 0
		self.tstartA = []
		self.tstopA = []
		self.expdataA = []
		self.ctsdataA = []
		self.ctssimA = []

		#results

		self.res = []
		return

	def loadDataAPAGILE(self, apfile):
		print('Loading data...' + apfile)
		self.apfile = apfile
		self.diml = 0
		with open(apfile, "r") as ins:
			for line in ins:
				if(line != ""):
					val = line.split(" ")
					expdata = float(val[2])
					if(expdata != 0):
						self.diml += 1

		self.tstartA = np.zeros(self.diml)
		self.tstopA = np.zeros(self.diml)
		self.expdataA = np.zeros(self.diml) #[cm2 s}
		self.ctsdataA = np.zeros(self.diml) #cts 

		nline = 0
		with open(apfile, "r") as ins:
			for line in ins:
				if(line != ""):
					val = line.split(" ")
					tstart  = float(val[0])
					tstop   = float(val[1])
					expdata = float(val[2])
					ctsdata = float(val[3]) #select the right column
					if(expdata != 0):
						self.tstartA[nline] = tstart
						self.tstopA[nline] = tstop
						self.expdataA[nline] = expdata
						self.ctsdataA[nline] = ctsdata
						nline += 1

							
	def calculateCtsFromConstantModel(self, ranal= 1, gasvalue=0, gal = 0.7, iso = 12., emin = 100., emax = 10000., gindex=2.1, source_theta=30, instrumentID = 0):
		if self.diml == 0:
			print('Error: data not loaded')
			return

		rate = EvalRates()
		bkg_ON, src_ON = rate.calculateRateWithoutExp(0, ranal, 0e-08,  gasvalue, gal, iso, emin, emax, gindex, source_theta, instrumentID)
		self.ctssimA = np.zeros(self.diml)
		ii = 0
		for exp in self.expdataA:
			ctsdatabkg = np.random.poisson(bkg_ON * exp)
			self.ctssimA[ii] = ctsdatabkg
			ii += 1
		return

	
	#load an ap4 file
	def loadnormalizedAP4(self, apfile):
		self.apfile = apfile
		n=0
		with open(apfile, "r") as ins:
			for line in ins:
				if(line != ""):
					val = line.split()
					if len(val) > 2:
						n = n + 1
		self.diml = n
			
		self.res = np.zeros((self.diml, self.ncols))
		self.tstartA = np.zeros(self.diml)
		self.tstopA = np.zeros(self.diml)
		self.expdataA = np.zeros(self.diml)
		self.ctsdataA = np.zeros(self.diml)
		
		n=0
		with open(apfile, "r") as ins:
			for line in ins:
				if(line != ""):
					if n == 0:
						n = 1
						continue
					val = line.split()
					if len(val) > 2:
						self.tstartA[n] = float(val[0])
						self.tstopA[n] = float(val[1])
						self.expdataA[n] = float(val[2])
						self.ctsdataA[n] = float(val[3])
						for i in range(0,self.ncols):
							self.res[n,i] = float(val[4+i])
						n = n + 1
			
		print("Loaded " + str(n) + " lines")


	#generation of AP4 file
	#Evaluation of the background.
	#Method 1: use gasvalue, gal and iso parameters. An analytical evaluation of the background is performed. This become the column rateWeightedMeanR4
	#Method 2: the mean background value is passed as parameter, using the rateMeanBkgExpected parameter. This become the default value if used
	#evalULalgorithm -> selection of the significance algorithm for the evaluation of the UL and sensitivity 
	# 1 -> Slima 2 -> Sa
	# if evalULalgorithm=1, ranalB=10
	def normalizeAP(self, apfile, ranal=2, rateMeanBkgExpected=-1, gasvalue=0.00054, gal=0.7, iso=10, emin=100, emax=10000, gindex=2.1, writevonmissesfiles=0, evalULalgorithm=1):
				
		#if self.diml == 0:
		self.loadDataAPAGILE(apfile)
		
		self.res = np.zeros((self.diml, self.ncols))
		
		mat = MathUtils()
		
		n=0
		for x in self.ctsdataA:
			n = n + 1
		
		nap = NormalizeAP()
		rate = EvalRates()
		
		rateBkgExpected = 0
		if rateMeanBkgExpected < 0:
			rate_bkg_ON, rate_src_ON = rate.calculateRateWithoutExp(verbose=0, ranalS=ranal, fluxsource=0e-08,  gasvalue=gasvalue, gal=gal, iso=iso, emin=emin, emax=emax, gindex=gindex, source_theta=30, instrumentID=0)
			rateBkgExpected = rate_bkg_ON
		else:
			rateBkgExpected = rateMeanBkgExpected
		
		rateWeightedMeanR1, rateWeightedMeanR2, rateWeightedMeanR3, rateWeightedMeanR4, rateWeightedMeanR1aa = nap.normalizeAB3(self.expdataA, self.ctsdataA, rateBkgExpected, self.res[:,0], self.res[:,1], self.res[:,2], self.res[:,3], self.res[:,4], self.res[:,5], self.res[:,6], self.res[:,7], self.res[:,8], self.res[:,9], self.res[:,10], self.res[:,11],self.res[:,12], self.res[:,13], self.res[:,14],self.res[:,15], self.res[:,16], self.res[:,46])
		
		fluxscalefactor=rate.getFluxScaleFactor(verbose=1, gindex=gindex, ranal=ranal, emin=emin, emax=emax)
		#/ 1.66
		n=0
		for e_i in self.expdataA:
			
			#flux -> flux_ratediffR4
			self.res[n,17] = self.res[n,13] / float(fluxscalefactor)
			
			#flux error -> flux_ratediffR4Error
			n_i = self.res[n,13] / float(fluxscalefactor) * e_i
			s_irms = mat.lowCountsError(n_i)
			s_i = s_irms / e_i
			self.res[n,18] = s_i 
			
			#Sa #formula paper Bulgarelli&Aboudan 
			if self.ctsdataA[n] > 0:
				aps = APSignificance()
				ctsB = rateBkgExpected * e_i
				#self.res[n,19] = -2 * np.log(np.exp(self.ctsdataA[n]-ctsBKG) * np.power(ctsBKG / self.ctsdataA[n], self.ctsdataA[n]))
				self.res[n,19] = aps.Sa(verbose=0, ctsB=ctsB, ctsTOT=self.ctsdataA[n])
			else:
				self.res[n,19] = 0
			
			#flux_rate
			self.res[n,20] = self.res[n,15] / float(fluxscalefactor)
			
			#flux_rateError
			n_i = self.res[n,15] / float(fluxscalefactor) * e_i
			s_irms = mat.lowCountsError(n_i)
			s_i = s_irms / e_i
			self.res[n,21] = s_i #flux error
			
			#calculation of expected background counts rate from BKG model
			#26:cts_rateWeightedMeanR4
			self.res[n,22] = rateBkgExpected * e_i
			
			#Slm #formula Li&Ma
			if self.ctsdataA[n] > 0:
				aps = APSignificance()
				N_off = rateBkgExpected * e_i
				N_on  = self.ctsdataA[n]
				# We are using the analytical model of the background. ranalB=10, that is the usual radius of analysis used for the evalutation of gal and iso coefficients with the MLE.
				lima = aps.lima(verbose=0, N_on = N_on, N_off = N_off, ranalS=ranal, ranalB=10)
				self.res[n,23] = lima
			else:
				self.res[n,23] = 0
				
			#da espandere qui, se si cambia il modo di valutare il bkg
			evaluatedBackground = self.res[n,22] #BKG da modello
			
			if evaluatedBackground > 0:
				ctsB = evaluatedBackground
				
				#1:LiMa
				#2: Sa
				algorithm=evalULalgorithm
				
				ctsSourceUL, SignUL = rate.calcCountsLimit(2, ctsB, ranal, algorithm=algorithm)
				rateUL = (ctsSourceUL / e_i)
				fluxUL = rateUL / fluxscalefactor
				self.res[n,24] = SignUL
				self.res[n,25] = ctsSourceUL
				self.res[n,26] = rateUL
				self.res[n,27] = fluxUL
				
				#sensitivity3
				ctsSourceSens3, SignSens3 = rate.calcCountsLimit(3, ctsB, ranal, algorithm=algorithm)
				rateSens3 = (ctsSourceSens3 / e_i)
				fluxSens3 = rateSens3 / fluxscalefactor
				self.res[n,28] = SignSens3
				self.res[n,29] = ctsSourceSens3
				self.res[n,30] = rateSens3
				self.res[n,31] = fluxSens3
				
				#sensitivity4
				ctsSourceSens4, SignSens4 = rate.calcCountsLimit(4, ctsB, ranal, algorithm=algorithm)
				rateSens4 = (ctsSourceSens4 / e_i)
				fluxSens4 = rateSens4 / fluxscalefactor
				self.res[n,32] = SignSens4
				self.res[n,33] = ctsSourceSens4
				self.res[n,34] = rateSens4
				self.res[n,35] = fluxSens4
				
				#sensitivity5
				ctsSourceSens5, SignSens5 = rate.calcCountsLimit(5, ctsB, ranal, algorithm=algorithm)
				rateSens5 = (ctsSourceSens5 / e_i)
				fluxSens5 = rateSens5 / fluxscalefactor
				self.res[n,36] = SignSens5
				self.res[n,37] = ctsSourceSens5
				self.res[n,38] = rateSens5
				self.res[n,39] = fluxSens5
				
				self.res[n, 40] = fluxscalefactor
				self.res[n, 41] = rateWeightedMeanR1
				self.res[n, 42] = rateWeightedMeanR2
				self.res[n, 43] = rateWeightedMeanR3
				self.res[n, 44] = rateWeightedMeanR4
				self.res[n, 45] = rateWeightedMeanR1aa

				self.res[n, 47] = self.tstopA[n] - self.tstartA[n]
			
			n = n + 1
			
		
		fileclean = open(apfile + ".ap4","w")

		
		print("* Res column number:  0:normrateAB11 1:normrateAB12 2:normrateAB13 3:normrateAB14 4:normrateAB21 5:normrateAB22 6:normrateAB23 7:normrateAB24 8:normrateAB11aa 9:normrateAB21aa 10:ratediffR1 11:ratediffR2 12:ratediffR3 13:ratediffR4 14:ratediffR1AA 15:rate 16:rateError 17:flux_ratediffR4 18:flux_ratediffR4Error 19:Sa 20:flux_rate 21:flux_rateError 22:cts_rateWeightedMeanR4 23:Slm 24:SignUL 25:ctsSourceUL 26:rateUL 27:fluxUL 28:SignSens3 29:ctsSourceSens3 30:rateSens3 31:fluxSens3 32:SignSens4 33:ctsSourceSens4 34:rateSens4 35:fluxSens4 36:SignSens5 36:ctsSourceSens5 38:rateSens5 39:fluxSens5 40:fluxscalefactor 41:rateWeightedMeanR1 42:rateWeightedMeanR2 43:rateWeightedMeanR3 44:rateWeightedMeanR4 45:rateWeightedMeanR1aa 46:expBasedRateError 47:temporalBinSize")
		
		header = "tstart tstop exp cts normrateAB11 normrateAB12 normrateAB13 normrateAB14 normrateAB21 normrateAB22 normrateAB23 normrateAB24 normrateAB11aa normrateAB21aa ratediffR1 ratediffR2 ratediffR3 ratediffR4 ratediffR1AA rate rateError flux_ratediffR4 flux_ratediffR4Error Sa flux_rate flux_rateError cts_rateWeightedMeanR4 Slm SignUL ctsSourceUL rateUL fluxUL SignSens3 ctsSourceSens3 rateSens3 fluxSens3 SignSens4 ctsSourceSens4 rateSens4 fluxSens4 SignSens5 ctsSourceSens5 rateSens5 fluxSens5 fluxscalefactor rateWeightedMeanR1 rateWeightedMeanR2 rateWeightedMeanR3 rateWeightedMeanR4 rateWeightedMeanR1aa expBasedRateError temporalBinSize"
		
		print("AP4 file column number: 0:tstart 1:tstop 2:exp[cm2s] 3:cts 4:normrateAB11 5:normrateAB12 6:normrateAB13 7:normrateAB14 8:normrateAB21 9:normrateAB22 10:normrateAB23 11:normrateAB24 12:normrateAB11aa 13:normrateAB21aa 14:ratediffR1 15:ratediffR2 16:ratediffR3 17:ratediffR4 18:ratediffR1AA 19:rate 20:rateError 21:flux_ratediffR4 22:flux_ratediffR4Error 23:Sa 24:flux_rate 25:flux_rateError 26:cts_rateWeightedMeanR4 27:Slm 28:SignUL 29:ctsSourceUL 30:rateUL 31:fluxUL 32:SignSens3 33:ctsSourceSens3 34:rateSens3 35:fluxSens3 36:SignSens4 37:ctsSourceSens4 38:rateSens4 39:fluxSens4 40:SignSens5 41:ctsSourceSens5 42:rateSens5 43:fluxSens5 44:fluxscalefactor 45:rateWeightedMeanR1 46:rateWeightedMeanR2 47:rateWeightedMeanR3 48:rateWeightedMeanR4 49:rateWeightedMeanR1aa 50:expBasedRateError 51:temporalBinSize")
		n = 0
		fileclean.write(header + "\n")
		for x in self.expdataA:
			line = str(self.tstartA[n]) + " " + str(self.tstopA[n]) + " " + str(self.expdataA[n]) + " " + str(int(self.ctsdataA[n]))
			for i in range(0,self.ncols):
				if i >= 0 and i <= 9:
					line += ' {:.2f}'.format(self.res[n,i])
				elif i > 9 and i <= 18:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i == 19:
					line += ' {:.2f}'.format(self.res[n,i])
				elif i >= 20 and i <= 21:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i >= 22 and i <= 25:
					line += ' {:.2f}'.format(self.res[n,i])
				elif i >= 26 and i <= 27:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i >= 28 and i <= 29:
					line += ' {:.2f}'.format(self.res[n,i])	
				elif i >= 30 and i <= 31:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i >= 32 and i <= 33:
					line += ' {:.2f}'.format(self.res[n,i])	
				elif i >= 34 and i <= 35:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i >= 36 and i <= 37:
					line += ' {:.2f}'.format(self.res[n,i])	
				elif i >= 38 and i <= 39:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i >= 40 and i <= 40:
					line += ' {:.2f}'.format(self.res[n,i])	
				elif i >= 41 and i <= 46:
					line += ' {:.2e}'.format(self.res[n,i])
				else:
					line += " " + str(self.res[n,i])
			line += "\n"
			fileclean.write(line)
			n = n + 1
		
		fileclean.close()

		print("* Write AP4 file: " + apfile + ".ap4")
		
		if(writevonmissesfiles == 1):
			self.writeVonMissesFile(apfile, 0)
			self.writeVonMissesFile(apfile, 1)
			self.writeVonMissesFile(apfile, 2)
			self.writeVonMissesFile(apfile, 3)
			self.writeVonMissesFile(apfile, 8)
			self.writeVonMissesFileRateA(apfile)
			self.writeVonMissesFileRateB(apfile)

		#write ap4r file
		filecleanr = open(apfile + ".ap4r","w")
		header = "tstart tstop exp cts rate rateError expBasedRateError flux_ratediffR4 flux_ratediffR4Error Sa Slm cts_rateWeightedMeanR4 rateWeightedMeanR4 ctsSourceUL rateUL fluxUL"
		n = 0
		filecleanr.write(header + "\n")
		for x in self.expdataA:
			line = str(self.tstartA[n]) + " " + str(self.tstopA[n]) + " " + str(self.expdataA[n]) + " " + str(int(self.ctsdataA[n]))
			#rate
			line += ' {:.2e}'.format(self.res[n,15])
			#rateError
			line += ' {:.2e}'.format(self.res[n,16])
			#expBasedRateError
			line += ' {:.2e}'.format(self.res[n,46])
			#flux_ratediffR4 
			line += ' {:.2e}'.format(self.res[n,17])
			#flux_ratediffR4Error
			line += ' {:.2e}'.format(self.res[n,18])
			#Sa
			line += ' {:.2f}'.format(self.res[n,19]) 
			#Slm 
			line += ' {:.2f}'.format(self.res[n,23])
			#cts_rateWeightedMeanR4 
			line += ' {:.2f}'.format(self.res[n,22])
			#rateWeightedMeanR4
			line += ' {:.2e}'.format(self.res[n,44])
			#ctsSourceUL
			line += ' {:.2f}'.format(self.res[n,25])
			#rateUL
			line += ' {:.2e}'.format(self.res[n,26])
			#fluxUL
			line += ' {:.2e}'.format(self.res[n,27])

			line += "\n"
			filecleanr.write(line)
			n = n + 1

		filecleanr.close()

		print('End normalisation')
		return
		

	def writeVonMissesFile(self, apfile, rescol):
		filevm = open(apfile + ".vm" + str(rescol),"w")
		n = 0
		for x in self.expdataA:
			if self.expdataA[n] > 0:
				line = str(self.tstartA[n] + (self.tstopA[n] - self.tstartA[n])/2.0) + " " + str(self.res[n,rescol]) + " " + str(np.sqrt(np.fabs(self.res[n,rescol]))) + "\n"
				filevm.write(line)
			n = n + 1

		filevm.close()
		print("* Write Von Misses file: tcenter rescol"+str(rescol)+" sqrt(abs(rescol"+str(rescol)+"))")

	"""Write Von Misses file rateA: tcenter rate rateError 
	"""
	def writeVonMissesFileRateA(self, apfile):
		filevm = open(apfile + ".vm16","w")
		rescol = 15
		n = 0
		for x in self.expdataA:
			if self.expdataA[n] > 0:
				line = str(self.tstartA[n] + (self.tstopA[n] - self.tstartA[n])/2.0) + " " + str(self.res[n,15]) + " " + str(self.res[n,16]) + "\n"
				filevm.write(line)
			n = n + 1

		filevm.close()
		print("* Write Von Misses file rateA (16): tcenter rate rateError")

	"""Write Von Misses file rateB: tcenter rate expBasedRateError
	"""
	def writeVonMissesFileRateB(self, apfile):
		filevm = open(apfile + ".vm46","w")
		rescol = 15
		n = 0
		for x in self.expdataA:
			if self.expdataA[n] > 0:
				line = str(self.tstartA[n] + (self.tstopA[n] - self.tstartA[n])/2.0) + " " + str(self.res[n,15]) + " " + str(self.res[n,46]) + "\n"
				filevm.write(line)
			n = n + 1
		
		filevm.close()
		print("* Write Von Misses file rateB (46): tcenter rate expBasedRateError")


	def fullAnalysis(self, apfilename, ranal=2, gasvalue=0.00054, analyzevm=-1, vonmissesthread=48, freqmin=1.0e-07, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800, vmnoise=1, gal=0.7, iso=10, emin=100, emax=10000, gindex=2.1, writevonmissesfiles=0, evalULalgorithm=1, rateMeanBkgExpected=-1):
		if analyzevm == 1:
			writevonmissesfiles = 1

		self.normalizeAP(apfile=apfilename, ranal=ranal, gasvalue=gasvalue, gal=gal, iso=iso, emin=emin, emax=emax, gindex=gindex, writevonmissesfiles=writevonmissesfiles, evalULalgorithm=evalULalgorithm, rateMeanBkgExpected=rateMeanBkgExpected)
		self.freqmin=float(freqmin)
		self.freqmax=float(freqmax)
		ls = MethodLS(self.tstartA, self.res, apfilename=apfilename)
		ls.scanLS(self.freqmin, self.freqmax, 0, 18)		

		if analyzevm == 1:
			vm = MethodVonMisses()
			vm.fullAnalysis(apfilename, vonmissesthread=vonmissesthread, freqmin=freqmin, freqmax=freqmax, vmnumax=vmnumax, ngridfreq=ngridfreq, tgridfreq=tgridfreq, vmnoise=vmnoise)

		print("End full analysis")


	def fullAnalysisLoadAP4(self, apfilename, analyzevm=-1, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800, vmnoise=1):
		self.loadnormalizedAP4(apfilename)
		self.freqmin=float(freqmin)
		self.freqmax=float(freqmax)
		ls = MethodLS(self.tstartA, self.res, apfilename=apfilename)
		ls.scanLS(self.freqmin, self.freqmax, 0, 18)
		
		if analyzevm == 1:
			vm = MethodVonMisses()
			vm.fullAnalysis(apfilename, vonmissesthread=vonmissesthread, freqmin=freqmin, freqmax=freqmax, vmnumax=vmnumax, ngridfreq=ngridfreq, tgridfreq=tgridfreq, vmnoise=vmnoise)

		print("End full analysis")
	
	def evaluateLS(self, apfile, rescol=3, plot=1, minfreq=0.5e-6, maxfreq=5e-6):
		self.apfile = apfile
		self.loadnormalizedAP4(apfile)
		ls = MethodLS(self.tstartA, self.res, apfilename=apfile)
		pls, pmax, maxf = ls.calculateLS(verbose=1, plot=plot, rescol=rescol, minfreq=minfreq, maxfreq=maxfreq)


	def evaluateLSexp(self, apfile, verbose=1, rescol=3, plot=1):
		self.apfile = apfile
		self.loadnormalizedAP4(apfile)
		ls = MethodLS(self.tstartA, self.res, apfilename=apfile)
		ls.calculateLSexp(verbose=verbose, plot=plot, rescol=rescol)

	def scanLS(self, apfile, verbose=1, plot=0):
		self.apfile = apfile
		self.loadnormalizedAP4(apfile)
		ls = MethodLS(self.tstartA, self.res, apfilename=apfile)
		ls.scanLS(self.freqmin, self.freqmax, 0, 18)
	
	def scanLSexp(self, apfile, verbose=1, plot=0):
		self.apfile = apfile
		self.loadnormalizedAP4(apfile)
		ls = MethodLS(self.tstartA, self.res, apfilename=apfile)
		for i in range(0, 18):
			ls.calculateLSexp(verbose=verbose, plot=plot, rescol=i)
