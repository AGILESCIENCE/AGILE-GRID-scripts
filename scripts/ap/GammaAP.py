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

import scipy.signal as signal
from astropy.stats import LombScargle
import string, os, sys
import numpy as np
import math
import matplotlib.pyplot as plt
from EvalRates import *
from MethodVonMisses import *
from APSignificance import *
from MathUtils import *

#gasvalue IGR = 0.00054
#gasvalue Vela = 0.00018
#TODO
#1) Includere la possibilita' di avere una LC di bkg indipendente da quella da analizzare, su cui determinare rate medio


class NormalizeAP:

	def __init__(self):
		return

	# Section 2
	#
	# [1] Robin Corbet, LAT Light Curve Analysis Aperture Photometry and Periodicity Searches, slides at COSPAR CBW 2010
	# [2] Gehrels, 1986, ApJ, 303, 336
	# [3] Kraft, Burrows, & Nousek, 1991, ApJ, 374, 344


	# Variance weighting corrected for low count statistic - Weighted Power Spectrum	 [1] [2]

	def normalizeAB3(self, expdataA, ctsdataA, rateBkgExpected, rateAB11, rateAB12, rateAB13, rateAB14, rateAB21, rateAB22, rateAB23, rateAB24, rateAB31, rateAB32, rateABR1, rateABR2, rateABR3, rateABR4, rateAAR5, rate, rate_error):
		sum1 = 0.0
		sum2 = 0.0
		sum3 = 0.0
		sum4 = 0.0
		sum5 = 0.0
		diml = len(expdataA)
		#rate = np.zeros(diml)
		#rate_error = np.zeros(diml)
		
		n=0
		for e_i in expdataA:
			
			n_i = ctsdataA[n] #cts
			
			rate_i = n_i / e_i #cts / cm2 s
			rate[n] = rate_i #cts / cm2 s
			
			s_ip = 0.5 + np.sqrt(n_i + 0.25)
			s_im = -0.5 + np.sqrt(n_i + 0.25)
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i
			rate_error[n] = s_i
			sum1 += rate_i / (s_i*s_i)
			sum2 += 1.0 / (s_i*s_i)
			sum3 += rate_i
			sum4 += n_i
			sum5 += e_i
			
			n = n + 1
			
		ratew_mean1 = sum1 / sum2 #cts / cm2 s
		ratew_mean2 = sum3 / n #cts / cm2 s
		ratew_mean3 = sum4 / sum5 #cts / cm2 s
		ratew_mean4 = rateBkgExpected #cts / cm2 s
		print("ratew_mean1:          %.3e"% ratew_mean1)
		print("ratew_mean2:          %.3e"% ratew_mean2)
		print("ratew_mean3:          %.3e"% ratew_mean3)
		print("ratew_mean4 bkgmodel: %.3e"% ratew_mean4)
		
		n=0
		sum1 = 0.0
		sum2 = 0.0
		for rate_i in rate:
			e_i = expdataA[n]
			ctspred_i = e_i * ratew_mean1 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			sum1 += rate_i / (sp_i*sp_i)
			sum2 += 1.0 / (sp_i*sp_i)
		
		ratew_mean1aa = sum1 / sum2 #cts / cm2 s
		print("ratew_mean1aa:        %.3e" %ratew_mean1aa)
		

		n=0
		for rate_i in rate:
			e_i = expdataA[n]
			s_i = rate_error[n]
			
			###############7.1.1
			#res column 0
			ctspred_i = e_i * ratew_mean1 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB11[n] = (rate_i - ratew_mean1) / (sp_i*sp_i)
			
			###############7.1.2
			#res column 1
			ctspred_i = e_i * ratew_mean2 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB12[n] = (rate_i - ratew_mean2) / (sp_i*sp_i)
			
			###############7.1.3
			#res column 2
			ctspred_i = e_i * ratew_mean3 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB13[n] = (rate_i - ratew_mean3) / (sp_i*sp_i)
			
			###############7.1.4
			#res column 3
			ctspred_i = e_i * ratew_mean4 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB14[n] = (rate_i - ratew_mean4) / (sp_i*sp_i)
			
			###############7.2.1
			#res column 4
			ctspred_i = e_i * ratew_mean1 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB21[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.2
			#res column 5
			ctspred_i = e_i * ratew_mean2 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB22[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.3
			#res column 6
			ctspred_i = e_i * ratew_mean3 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB23[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.4
			#res column 7
			ctspred_i = e_i * ratew_mean4 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB24[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.3.1
			#rateAB31[n] = (rate_i - ratew_mean1) / (s_i*s_i)
			
			###############7.3.2
			#rateAB32[n] = (rate_i - ratew_mean2) / (s_i*s_i)
			
			###############7.3.3
			#rateAB33[n] = (rate_i - ratew_mean3) / (s_i*s_i)
			
			###############7.3.4
			#rateAB34[n] = (rate_i - ratew_mean4) / (s_i*s_i)
			
			###############7.1.1aa
			#res column 8
			ctspred_i = e_i * ratew_mean1aa #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB31[n] = (rate_i - ratew_mean1aa) / (sp_i*sp_i)
			
			###############7.2.1aa
			#res column 9
			ctspred_i = e_i * ratew_mean1aa #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB32[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.4.1
			#rateAB41[n] = (rate_i ) / (s_i*s_i)
			
			###############rate
			#res column 10
			rateABR1[n] = (rate_i - ratew_mean1) 
			
			###############rate
			#res column 11
			rateABR2[n] = (rate_i - ratew_mean2)
			
			###############rate
			#res column 12
			rateABR3[n] = (rate_i - ratew_mean3)
			
			###############rate
			#res column 13
			rateABR4[n] = (rate_i - ratew_mean4)
			
			#res column 14
			rateAAR5[n] = (rate_i - ratew_mean1aa)
			
			n = n + 1
		
		return


#########################################
class GammaAP:

	def __init__(self):
		self.ncols = 38
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

	
	#load an ap3 file
	def loadnormalizedAP3(self, ap3file):
		self.apfile = ap3file
		n=0
		with open(ap3file, "r") as ins:
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
		with open(ap3file, "r") as ins:
			for line in ins:
				if(line != ""):
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


	#generation of AP3 file
	def normalizeAP(self, apfile, ranal=2, gasvalue=0.00054, gal=0.7, iso=10, emin=100, emax=10000, gindex=2.1, writevonmissesfiles=0):
				
		#if self.diml == 0:
		self.loadDataAPAGILE(apfile)
		
		self.res = np.zeros((self.diml, self.ncols))
		
		mat = MathUtils()
		
		n=0
		for x in self.ctsdataA:
			n = n + 1
		
		nap = NormalizeAP()
		
		rate = EvalRates()
		rate_bkg_ON, rate_src_ON = rate.calculateRateWithoutExp(verbose=0, ranalS=ranal, fluxsource=0e-08,  gasvalue=gasvalue, gal=gal, iso=iso, emin=emin, emax=emax, gindex=gindex, source_theta=30, instrumentID=0)
		rateBkgExpected = rate_bkg_ON
		
		nap.normalizeAB3(self.expdataA, self.ctsdataA, rateBkgExpected, self.res[:,0], self.res[:,1], self.res[:,2], self.res[:,3], self.res[:,4], self.res[:,5], self.res[:,6], self.res[:,7], self.res[:,8], self.res[:,9], self.res[:,10], self.res[:,11],self.res[:,12], self.res[:,13], self.res[:,14],self.res[:,15], self.res[:,16])
		
		fluxscalefactor=rate.getFluxScaleFactor(verbose=1, gindex=gindex, ranal=ranal, emin=emin, emax=emax)
		#/ 1.66
		n=0
		for e_i in self.expdataA:
			
			#flux -> flux_ratediffR4
			self.res[n,17] = self.res[n,13] / float(fluxscalefactor)
			
			#flux error -> flux_ratediffR4_error
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
				self.res[n,19] = -1
			
			#flux rate
			self.res[n,20] = self.res[n,15] / float(fluxscalefactor)
			
			#flux rate error
			n_i = self.res[n,15] / float(fluxscalefactor) * e_i
			s_irms = mat.lowCountsError(n_i)
			s_i = s_irms / e_i
			self.res[n,21] = s_i #flux error
			
			#calculation of expected background counts rate
			#26:cts_expBKG4
			self.res[n,22] = rateBkgExpected * e_i
			
			#Slm #formula Li&Ma
			if self.ctsdataA[n] > 0:
				aps = APSignificance()
				N_off = rateBkgExpected * e_i
				N_on  = self.ctsdataA[n]
				lima = aps.lima(verbose=0, N_on = N_on, N_off = N_off, ranalS=ranal)
				self.res[n,23] = lima
			else:
				self.res[n,23] = -1
				
			if self.ctsdataA[n] > 0:
				ctsB = self.res[n,22]
				
				#1:LiMa
				#2: Sa
				algorithm=1
				
				N_sourceUL, SignUL = rate.calcCountsLimit(2, ctsB, ranal, algorithm=algorithm)
				rateUL = (N_sourceUL / e_i)
				fluxUL = rateUL / fluxscalefactor
				self.res[n,24] = SignUL
				self.res[n,25] = N_sourceUL
				self.res[n,26] = rateUL
				self.res[n,27] = fluxUL
				
				#sensitivity4
				N_sourceSens4, SignSens4 = rate.calcCountsLimit(4, ctsB, ranal, algorithm=algorithm)
				rateSens4 = (N_sourceSens4 / e_i)
				fluxSens4 = rateSens4 / fluxscalefactor
				self.res[n,28] = SignSens4
				self.res[n,29] = N_sourceSens4
				self.res[n,30] = rateSens4
				self.res[n,31] = fluxSens4
				
				#sensitivity5
				N_sourceSens5, SignSens5 = rate.calcCountsLimit(5, ctsB, ranal, algorithm=algorithm)
				rateSens5 = (N_sourceSens5 / e_i)
				fluxSens5 = rateSens5 / fluxscalefactor
				self.res[n,32] = SignSens5
				self.res[n,33] = N_sourceSens5
				self.res[n,34] = rateSens5
				self.res[n,35] = fluxSens5
				
				self.res[n, 36] = fluxscalefactor
				self.res[n, 37] = rateBkgExpected
			
			n = n + 1
			
		
		fileclean = open(apfile + ".ap3","w")

		
		print("* AP3 file res column number:   tstart tstop exp[cm2s] cts 0:normAB11 1:normAB12 2:normAB13 3:normAB14 4:normAB21 5:normAB22 6:normAB23 7:normAB24 8:normAB11aa 9:normAB21aa 10:ratediffR1 11:ratediffR2 12:ratediffR3 13:ratediffR4 14:ratediffR1AA 15:rate 16:rate_error 17:flux_ratediffR4 18:flux_ratediffR4_error 19:Sa 20:flux_rate 21:flux_rate_error 22:cts_expBKG4 23:Slm 24:SignUL 25:N_sourceUL 26:rateUL 27:fluxUL 28:SignSens4 29:N_sourceSens4 30:rateSens4 31:fluxSens4 32:SignSens5 33:N_sourceSens5 34:rateSens5 35:fluxSens5 36:fluxscalefactor")
		header = "tstart tstop exp cts normAB11 normAB12 normAB13 normAB14 normAB21 normAB22 normAB23 normAB24 normAB11aa normAB21aa ratediffR1 ratediffR2 ratediffR3 ratediffR4 ratediffR1AA rate rate_error flux_ratediffR4 flux_ratediffR4_error Sa flux_rate flux_rate_error cts_expBKG4 Slm SignUL N_sourceUL rateUL fluxUL SignSens4 N_sourceSens4 rateSens4 fluxSens4 SignSens5 N_sourceSens5 rateSens5 fluxSens5 fluxscalefactor"
		
		print("AP3 file column number: 0:tstart 1:tstop 2:exp[cm2s] 3:cts 4:normAB11 5:normAB12 6:normAB13 7:normAB14 8:normAB21 9:normAB22 10:normAB23 11:normAB24 12:normAB11aa 13:normAB21aa 14:ratediffR1 15:ratediffR2 16:ratediffR3 17:ratediffR4 18:ratediffR1AA 19:rate 20:rate_error 21:flux_ratediffR4 22:flux_ratediffR4_error 23:Sa 24:flux_rate 25:flux_rate_error 26:cts_expBKG4 27:Slm 28:SignUL 29:N_sourceUL 30:rateUL 31:fluxUL 32:SignSens4 33:N_sourceSens4 34:rateSens4 35:fluxSens4 36:SignSens5 37:N_sourceSens5 38:rateSens5 39:fluxSens5 40:fluxscalefactor")
		
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
				elif i > 19 and i <= 21:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i > 21 and i <= 25:
					line += ' {:.2f}'.format(self.res[n,i])
				elif i > 25 and i <= 27:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i > 27 and i <= 29:
					line += ' {:.2f}'.format(self.res[n,i])	
				elif i > 29 and i <= 31:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i > 31 and i <= 33:
					line += ' {:.2f}'.format(self.res[n,i])	
				elif i > 33 and i <= 35:
					line += ' {:.2e}'.format(self.res[n,i])
				elif i == 36:
					line += ' {:.4f}'.format(self.res[n,i])
				elif i == 37:
					line += ' {:.2e}'.format(self.res[n,i])
				else:
					line += " " + str(self.res[n,i])
			line += "\n"
			fileclean.write(line)
			n = n + 1
		
		fileclean.close()
		
		if(writevonmissesfiles == 1):
			self.writeVonMissesFile(apfile, 0)
			self.writeVonMissesFile(apfile, 1)
			self.writeVonMissesFile(apfile, 2)
			self.writeVonMissesFile(apfile, 3)
			self.writeVonMissesFile(apfile, 8)
		print('End normalisation')
		return
		

	def writeVonMissesFile(self, apfile, ii):
		filevm = open(apfile + ".vm" + str(ii),"w")
		n = 0
		for x in self.expdataA:
			line = str(self.tstartA[n] + (self.tstopA[n] - self.tstartA[n])/2.0) + " " + str(self.res[n,ii]) + " " + str(np.sqrt(np.fabs(self.res[n,ii]))) + "\n"
			filevm.write(line)
			n = n + 1

		filevm.close()
		print("* Write Von Misses file: tcenter  "+str(ii)+":rate 1:rate_var")


	def calculateLS(self, verbose=0, plot=1, rescol=0, minfreq=-1, maxfreq=-1):
		#normalization standard model log psd
		#, normalization='standard'
		ls = LombScargle(self.tstartA, self.res[:,rescol])
		if minfreq != -1:
			#10e-6, 10e-4 , samples_per_peak=10
			frequency, power = ls.autopower(minimum_frequency=minfreq, maximum_frequency=maxfreq)
		else:
			frequency, power = ls.autopower()
		pls = ls.false_alarm_probability(power.max(), method='baluev')
		ind = power.argmax()
		maxf = frequency[ind]
		daymax = 1 / maxf / 86400.

		if verbose == 1:
			print('Max value col %2d: FAR %.5e / power %.4f / frequency %.5e Hz (%.5f days)' %(rescol, pls, power.max(), maxf, daymax) )

		if plot > 0:
			plt.subplot(1, 1, 1)
			plt.plot(frequency, power)
			plt.gca().set_xscale("log")
			plt.gca().set_yscale("log")
			#plt.subplot(3, 1, 2)
			#best_frequency = frequency[np.argmax(power)]
			#t_fit = np.linspace(0, 1)
			#y_fit = LombScargle(tstartA, ctsdataA).model(t_fit, best_frequency)
			#plt.plot(t_fit, y_fit)
		if plot == 1:
			plt.show()
		if plot == 2:
			plt.savefig(self.apfile + ".apres.ls."+str(rescol)+".png", format='png', dpi=1200)

		return pls, power.max(), maxf

	def scanLS(self, minfreq=-1, maxfreq=-1, rangemin=0, rangemax=23):
		print("freqmin %.5e"%minfreq)
		print("freqmax %.5e"%maxfreq)
		fileclean = open(self.apfile + ".apres","w")
		for i in range(rangemin,rangemax):
			pls, pmax, maxf = self.calculateLS(0 ,0, i, minfreq, maxfreq)
			daymax = 1 / maxf / 86400.
			#print("LS " + str(i) + ' ' + str(pls) + ' ' + str(pmax) + ' ' + str(maxf) + ' ' + str(daymax))
			print('Max value %2d: FAR %.5e / power %.4f / frequency %.5e Hz (%.5f days)' %(i, pls, pmax, maxf, daymax) )
			fileclean.write("LS " + str(i) + ' ' + str(pls) + ' ' + str(pmax) + ' ' + str(maxf) + ' ' + str(daymax) + "\n")
		fileclean.close()
		return

	def calculateLSexp(self, verbose=0, plot=1):
		ls = LombScargle(self.tstartA, self.expdataA)
		frequency, power = ls.autopower()
		pls = ls.false_alarm_probability(power.max(), method='baluev')
		ind = power.argmax()
		maxf = frequency[ind]
		daymax = 1 / maxf / 86400.
		if verbose == 1:
			print('pls ' + str(pls) + ' with power of ' + str(power.max()) + ' at frequency ' + str(maxf) + ' Hz (' + str(daymax) + ') days' )

		if plot == 1:
			plt.subplot(1, 1, 1)
			plt.plot(frequency, power)
			plt.gca().set_xscale("log")
			plt.gca().set_yscale("log")
			#plt.subplot(3, 1, 2)
			#best_frequency = frequency[np.argmax(power)]
			#t_fit = np.linspace(0, 1)
			#y_fit = LombScargle(tstartA, ctsdataA).model(t_fit, best_frequency)
			#plt.plot(t_fit, y_fit)
			plt.show()

	def plotLS(self, ap3file, i):
		self.apfile = ap3file
		self.loadnormalizedAP3(ap3file)
		pls, pmax, maxf = self.calculateLS(1, 2, int(i), 0.5e-6, 5e-6)


	def fullAnalysis(self, apfilename, ranal=2, gasvalue=0.00054, analyzevm=-1, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800):
		self.normalizeAP(apfilename, ranal, gasvalue)
		self.freqmin=float(freqmin)
		self.freqmax=float(freqmax)
		
		self.scanLS(self.freqmin, self.freqmax, 0, 18)
		#self.scanLS(-1, -1, 0, 19)

		if analyzevm == 1:
			vm = MethodVonMisses()
			vm.fullAnalysis(apfilename, vonmissesthread, freqmin, freqmax, vmnumax, ngridfreq, tgridfreq)


	def fullAnalysisLoadAP3(self, ap3filename, analyzevm=-1, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800):
		self.loadnormalizedAP3(ap3filename)
		self.freqmin=float(freqmin)
		self.freqmax=float(freqmax)
		self.scanLS(self.freqmin, self.freqmax, 0, 18)
		
		if analyzevm == 1:
			vm = MethodVonMisses()
			vm.fullAnalysis(ap3filename, vonmissesthread, freqmin, freqmax, vmnumax, ngridfreq, tgridfreq)
