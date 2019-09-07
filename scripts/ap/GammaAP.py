import scipy.signal as signal
from astropy.stats import LombScargle
import string, os, sys
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from EvalRates import *

#gasvalue IGR = 0.00054
#gasvalue Vela = 0.00018
#TODO
#1) Includere la possibilità di avere una LC di bkg indipendente da quella da analizzare, su cui determinare rate medio


class NormalizeAP:

	def __init__(self):
		return

	# https://docs.google.com/document/d/1FLSHiYUU1Lcl2JFkwAMKMU0GG9GbX4Ou0jUADRYIPr8/edit?usp=sharing
	# Section 2
	#
	# [1] Robin Corbet, LAT Light Curve Analysis: Aperture Photometry and Periodicity Searches, slides at COSPAR CBW 2010
	# [2] Gehrels, 1986, ApJ, 303, 336
	# [3] Kraft, Burrows, & Nousek, 1991, ApJ, 374, 344


	# Variance weighting corrected for low count statistic - Weighted Power Spectrum	 [1] [2]
	#se ctssimA e' passato, calcola il rate medio sul dato simulato
	def normalizeDef(self, expdataA, ctsdataA, ratew, rate, s, ctssimA):
		sum1 = 0.0
		sum2 = 0.0
		diml = len(expdataA)

		n=0
		for e_i in expdataA:
			if len(ctssimA) > 0:
				nb_i = ctssimA[n]
			else:
				nb_i = ctsdataA[n]
			n_i = ctsdataA[n]

			#calcola il rate medio sul dato vero o sul dato simulato
			rate_i = nb_i / e_i
			s_ip = 0.5 + np.sqrt(nb_i + 0.25)
			s_im = -0.5 + np.sqrt(nb_i + 0.25)
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i
			sum1 += rate_i / (s_i*s_i)
			sum2 += 1.0 / (s_i*s_i)

			#calcola rate e errore su dato vero
			rate_i = n_i / e_i
			rate[n] = rate_i
			s_ip = 0.5 + np.sqrt(n_i + 0.25)
			s_im = -0.5 + np.sqrt(n_i + 0.25)
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i
			s[n] = s_i

			n = n + 1

		#calcolo del rate medio
		ratew_mean = sum1 / sum2

		n=0
		for rate_i in rate:
			s_i =  s[n]
			if s_i != 0.0:
				#def_ratew << (rate_i - ratew_mean) / s_i
				ratew[n] = (rate_i - ratew_mean) / (s_i*s_i)
			else:
				ratew[n] = 0.0

			n = n + 1

		return

	def normalizeAA(self, expdataA, ctsdataA, rateAA, ctssimA):
		#AA
		sum1 = 0.0
		sum2 = 0.0
		sum3 = 0.0
		diml = len(expdataA)
		rate = np.zeros(diml)

		n=0
		for e_i in expdataA:
			if len(ctssimA) > 0:
				nb_i = ctssimA[n]
			else:
				nb_i = ctsdataA[n]
			n_i = ctsdataA[n]

			rate_i = nb_i / e_i
			#method 1
			s_ip = 0.5 + np.sqrt(nb_i + 0.25)
			s_im = -0.5 + np.sqrt(nb_i + 0.25)
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i

			sum1 += rate_i / (s_i*s_i)
			sum2 += 1.0 / (s_i*s_i)
			sum3 += rate_i

			rate_i = n_i / e_i
			rate[n] = rate_i

			n=n+1

		rate_mean = sum3 / n
		ratew_mean = sum1 / sum2

		n=0
		for rate_i in rate:
			ratepred_i = expdataA[n] * rate_mean
			s_i = np.sqrt(ratepred_i) / expdataA[n]
			#aa_ratew << (rate_i - ratew_mean) / s_i
			if s_i != 0.0:
				rateAA[n] = (rate_i - ratew_mean) / (s_i*s_i)
			else:
				rateAA[n] = 0.0

			n = n + 1

		return

	def normalizeAB2(self, expdataA, ctsdataA, rateAB1, rateAB2, rateAB3, ctssimA):
		sum1 = 0.0
		sum2 = 0.0
		sum3 = 0.0
		sum4 = 0.0
		sum5 = 0.0
		diml = len(expdataA)
		rate = np.zeros(diml)

		ab1_ratew = []
		ab2_ratew = []
		ab3_ratew = []
		n=0
		for e_i in expdataA:
			if len(ctssimA) > 0:
				nb_i = ctssimA[n]
			else:
				nb_i = ctsdataA[n]
			n_i = ctsdataA[n]

			rate_i = nb_i / e_i
			s_ip = 0.5 + np.sqrt(nb_i + 0.25)
			s_im = -0.5 + np.sqrt(nb_i + 0.25)
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i
			sum1 += rate_i / (s_i*s_i)
			sum2 += 1.0 / (s_i*s_i)
			sum3 += rate_i
			sum4 += nb_i
			sum5 += e_i

			rate_i = n_i / e_i
			rate[n] = rate_i
			s_ip = 0.5 + np.sqrt(n_i + 0.25)
			s_im = -0.5 + np.sqrt(n_i + 0.25)
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i
			n = n + 1

		ratew_mean1 = sum1 / sum2
		ratew_mean2 = sum3 / n
		ratew_mean3 = sum4 / sum5
		n=0
		for rate_i in rate:
			e_i = expdataA[n]

			ratepred_i = e_i * ratew_mean1
			s_i = np.sqrt(ratepred_i) / e_i
			#ab1_ratew << (rate_i - ratew_mean1) / s_i
			rateAB1[n] = (rate_i - ratew_mean1) / (s_i*s_i)

			ratepred_i = e_i * ratew_mean2
			s_i = np.sqrt(ratepred_i) / e_i
			#ab2_ratew << (rate_i - ratew_mean2) / s_i
			rateAB2[n] = (rate_i - ratew_mean2) / (s_i*s_i)

			#red text method
			ratepred_i = e_i * ratew_mean3
			s_i = np.sqrt(ratepred_i) / e_i
			#ab3_ratew << (rate_i - ratew_mean3) / s_i
			rateAB3[n] = (rate_i - ratew_mean3) / (s_i*s_i)

			n = n + 1

		return


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
			#0
			ctspred_i = e_i * ratew_mean1 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB11[n] = (rate_i - ratew_mean1) / (sp_i*sp_i)
			
			###############7.1.2
			#1
			ctspred_i = e_i * ratew_mean2 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB12[n] = (rate_i - ratew_mean2) / (sp_i*sp_i)
			
			###############7.1.3
			#2
			ctspred_i = e_i * ratew_mean3 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB13[n] = (rate_i - ratew_mean3) / (sp_i*sp_i)
			
			###############7.1.4
			#3
			ctspred_i = e_i * ratew_mean4 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB14[n] = (rate_i - ratew_mean4) / (sp_i*sp_i)
			
			###############7.2.1
			#4
			ctspred_i = e_i * ratew_mean1 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB21[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.2
			#5
			ctspred_i = e_i * ratew_mean2 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB22[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.3
			#6
			ctspred_i = e_i * ratew_mean3 #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB23[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.2.4
			#7
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
			#8
			ctspred_i = e_i * ratew_mean1aa #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB31[n] = (rate_i - ratew_mean1aa) / (sp_i*sp_i)
			
			###############7.2.1aa
			#9
			ctspred_i = e_i * ratew_mean1aa #cts
			sp_i = np.sqrt(ctspred_i) / e_i
			rateAB32[n] = (rate_i) / (sp_i*sp_i)
			
			###############7.4.1
			#rateAB41[n] = (rate_i ) / (s_i*s_i)
			
			###############rate
			#10
			rateABR1[n] = (rate_i - ratew_mean1) 
			
			###############rate
			#11
			rateABR2[n] = (rate_i - ratew_mean2)
			
			###############rate
			#«12
			rateABR3[n] = (rate_i - ratew_mean3)
			
			###############rate
			#13
			rateABR4[n] = (rate_i - ratew_mean4)
			
			#14
			rateAAR5[n] = (rate_i - ratew_mean1aa)
			
			n = n + 1
		
		return


#########################################
class GammaAP:

	def __init__(self):
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
		self.expdataA = np.zeros(self.diml)
		self.ctsdataA = np.zeros(self.diml)

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

	def calculateCtsFromConstantModel(self, ranal= 1, gasvalue=0, gal = 0.7, iso = 12., emin = 100., emax = 10000., instrumentID = 0):
		if self.diml == 0:
			print('Error: data not loaded')
			return

		rate = EvalRates()
		bkg_ON, src_ON = rate.calculateRateWithoutExp(0, ranal, 0e-08,  gasvalue, gal, iso, emin, emax, instrumentID)
		self.ctssimA = np.zeros(self.diml)
		ii = 0
		for exp in self.expdataA:
			ctsdatabkg = np.random.poisson(bkg_ON * exp)
			self.ctssimA[ii] = ctsdatabkg
			ii += 1
		return

	#load an ap2 file
	def loadnormalizedAP2(self, ap2file):
		self.apfile = ap2file
		n=0
		with open(ap2file, "r") as ins:
			for line in ins:
				if(line != ""):
					val = line.split()
					if len(val) > 2:
						n = n + 1
		self.diml = n

		self.res = np.zeros((self.diml, 14))
		self.tstartA = np.zeros(self.diml)
		self.tstopA = np.zeros(self.diml)
		self.expdataA = np.zeros(self.diml)
		self.ctsdataA = np.zeros(self.diml)
		n=0
		with open(ap2file, "r") as ins:
			for line in ins:
				if(line != ""):
					val = line.split()
					if len(val) > 2:
							self.tstartA[n] = float(val[0])
							self.tstopA[n] = float(val[1])
							self.expdataA[n] = float(val[2])
							self.ctsdataA[n] = float(val[3])
							for i in range(0,12):
								self.res[n,i] = float(val[4+i])
							n = n + 1

		n=0
		for x in self.ctsdataA:
			self.res[n,12] = self.ctsdataA[n]
			self.res[n,13] = self.expdataA[n]
			n = n + 1

		print("Loaded " + str(n) + " lines")

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
			
		self.res = np.zeros((self.diml, 23))
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
						for i in range(0,23):
							self.res[n,i] = float(val[4+i])
						n = n + 1
			
		print("Loaded " + str(n) + " lines")

	def normalizeAP2(self, apfile):

		if self.diml == 0:
			self.loadDataAPAGILE(apfile)
			#return

		self.res = np.zeros((self.diml, 14))
		n=0
		for x in self.ctsdataA:
			self.res[n,12] = self.ctsdataA[n]
			self.res[n,13] = self.expdataA[n]
			n = n + 1

		nap = NormalizeAP()

		#0 rate cts/exp
		#1 variance of the rate
		#2 def_ratew
		nap.normalizeDef(self.expdataA, self.ctsdataA, self.res[:,2], self.res[:,0], self.res[:,1], [])

		#3 aa_ratew
		nap.normalizeAA(self.expdataA, self.ctsdataA, self.res[:,3], [])

		#4 ab1_ratew
		#5 ab2_ratew
		#6 ab3_ratew
		nap.normalizeAB2(self.expdataA, self.ctsdataA, self.res[:,4], self.res[:,5], self.res[:,6], [])

		if len(self.ctssimA) == 0:
			self.calculateCtsFromConstantModel(2, 0.00054, 1, 15, 100, 10000, 0)

		if len(self.ctssimA) > 0:
			#0 rate cts/exp
			#1 variance of the rate
			#7 def_ratewB
			nap.normalizeDef(self.expdataA, self.ctsdataA, self.res[:,7], self.res[:,0], self.res[:,1], self.ctssimA)
			#8 aa_ratewB
			nap.normalizeAA(self.expdataA, self.ctsdataA, self.res[:,8], self.ctssimA)
			#9 ab1_ratewB
			#10 ab2_ratewB
			#11 ab3_ratewB
			nap.normalizeAB2(self.expdataA, self.ctsdataA, self.res[:,9], self.res[:,10], self.res[:,11], self.ctssimA)

		fileclean = open(apfile + ".ap2","w")

		n = 0
		for x in self.expdataA:
			line = str(self.tstartA[n]) + " " + str(self.tstopA[n]) + " " + str(self.expdataA[n]) + " " + str(int(self.ctsdataA[n]))
			for i in range(0,12):
				line += " " + str(self.res[n,i])
			line += "\n"
			fileclean.write(line)
			n = n + 1

		fileclean.close()
		print("* Write ap2 file: tstart tstop exp[cm2 s] cts 0:rate[cts/exp] 1:rate_var 2:def_ratew 3:aa_ratew 4:ab1_ratew 5:ab2_ratew 6:ab3_ratew 7:def_ratewB 8:aa_ratewB 9:ab1_ratewB 10:ab2_ratewB 11:ab3_ratewB")
		print("  B means normalization for a background model: calculateCtsFromConstantModel(2, 0.00054, 1, 30, 100, 10000, 0)")

		self.writeVonMisses(apfile, 2)
		self.writeVonMisses(apfile, 3)
		self.writeVonMisses(apfile, 5)
		self.writeVonMisses(apfile, 7)
		self.writeVonMisses(apfile, 10)
		print('End normalisation')
		return


	def normalizeAP3(self, apfile, ranal=2, gasvalue=0.00054, gal=0.7, iso=10, emin=100, emax=10000):
				
		#if self.diml == 0:
		self.loadDataAPAGILE(apfile)
		
		self.res = np.zeros((self.diml, 23))
		
		n=0
		for x in self.ctsdataA:
			n = n + 1
		
		nap = NormalizeAP()
		
		rate = EvalRates()
		bkg_ON, src_ON = rate.calculateRateWithoutExp(0, ranal, 0e-08,  gasvalue, gal, iso, emin, emax, 0)
		rateBkgExpected = bkg_ON
		#print(bkg_ON)
		
		nap.normalizeAB3(self.expdataA, self.ctsdataA, rateBkgExpected, self.res[:,0], self.res[:,1], self.res[:,2], self.res[:,3], self.res[:,4], self.res[:,5], self.res[:,6], self.res[:,7], self.res[:,8], self.res[:,9], self.res[:,10], self.res[:,11],self.res[:,12], self.res[:,13], self.res[:,14],self.res[:,15], self.res[:,16])
		
		
		fluxscalefactor=rate.getFluxScaleFactor(0, ranal, emin, emax)
		n=0
		for e_i in self.expdataA:
			
			#flux
			self.res[n,17] = self.res[n,13] / float(fluxscalefactor)
			
			#flux error
			n_i = self.res[n,13] / float(fluxscalefactor) * e_i
			if n_i + 0.25 > 0:
				s_ip = 0.5 + np.sqrt(n_i + 0.25)
				s_im = -0.5 + np.sqrt(n_i + 0.25)
			else:
				s_ip = 0
				s_im = 0
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i
			self.res[n,18] = s_i 
			
			#TS
			if self.ctsdataA[n] > 0:
				self.res[n,19] = -2 * np.log(np.exp(self.ctsdataA[n]-rateBkgExpected * e_i) * np.power(rateBkgExpected * e_i / self.ctsdataA[n], self.ctsdataA[n]))
			else:
				self.res[n,19] = -1
			
			#flux rate
			self.res[n,20] = self.res[n,15] / float(fluxscalefactor)
			
			#flux rate error
			n_i = self.res[n,15] / float(fluxscalefactor) * e_i
			s_ip = 0.5 + np.sqrt(n_i + 0.25)
			s_im = -0.5 + np.sqrt(n_i + 0.25)
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			s_i = s_irms / e_i
			self.res[n,21] = s_i #flux error
			
			#calculation of expected background counts rate
			#26:cts_expBKG4
			self.res[n,22] = rateBkgExpected * e_i
			
			n = n + 1
			
		
		fileclean = open(apfile + ".ap3","w")

		n = 0
		for x in self.expdataA:
			line = str(self.tstartA[n]) + " " + str(self.tstopA[n]) + " " + str(self.expdataA[n]) + " " + str(int(self.ctsdataA[n]))
			for i in range(0,23):
				line += " " + str(self.res[n,i])
			line += "\n"
			fileclean.write(line)
			n = n + 1
		
		fileclean.close()
		print("* Write ap3 file: tstart tstop exp[cm2 s] cts 0:normAB11 1:normAB12 2:normAB13 3:normAB14 4:normAB21 5:normAB22 6:normAB23 7:normAB24 8:normAB11aa 9:normAB21aa 10:ratediffR1 11:ratediffR2 12:ratediffR3 13:ratediffR4 14:ratediffR1AA 15:rate 16:rate_error 17:flux_ratediffR4 18:flux_ratediffR4_error 19:TS  20:flux_rate 21:flux_rate_error 22:cts_expBKG4")
		print("AP3 file column numbers: 0:tstart 1:tstop 2:exp[cm2 s] 3:cts 4:normAB11 5:normAB12 6:normAB13 7:normAB14 8:normAB21 9:normAB22 10:normAB23 11:normAB24 12:normAB11aa 13:normAB21aa 14:ratediffR1 15:ratediffR2 16:ratediffR3 17:ratediffR4 18:ratediffR1AA 19:rate 20:rate_error 21:flux_ratediffR4 22:flux_ratediffR4_error 23:TS 24:flux_rate 25:flux_rate_error 26:cts_expBKG4")
		
		self.writeVonMisses(apfile, 0)
		self.writeVonMisses(apfile, 1)
		self.writeVonMisses(apfile, 2)
		self.writeVonMisses(apfile, 3)
		self.writeVonMisses(apfile, 8)
		print('End normalisation')
		return
		

	def writeVonMisses(self, apfile, ii):
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
		#rescol
		#0 rate cts/exp -> NO
		#1 variance of the rate -> NO
		#2 def_ratew -> NO
		#3 aa_ratew
		#4 ab1_ratew
		#5 ab2_ratew
		#6 ab3_ratew
		#7 def_ratewB
		#8 aa_ratewB
		#9 ab1_ratewB
		#10 ab2_ratewB
		#11 ab3_ratewB
		#12 cts
		#13 exp
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

	def scanLS(self, minfreq=-1, maxfreq=-1, rangemin=0, rangemax=14):
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

	def plotLS(self, ap2file, i):
		self.apfile = ap2file
		self.loadnormalizedAP2(ap2file)
		pls, pmax, maxf = self.calculateLS(1, 2, int(i), 0.5e-6, 5e-6)

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

	def fullAnalysis2(self, apfilename, analyzevm=-1, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800):
		self.normalizeAP2(apfilename)
		self.freqmin=float(freqmin)
		self.freqmax=float(freqmax)
		self.vmnumax=float(vmnumax)

		self.scanLS(self.freqmin, self.freqmax, 0, 14)

		if analyzevm == 1:
			self.runVomMisses(vonmissesthread, 2)
			self.runVomMisses(vonmissesthread, 3)
			self.runVomMisses(vonmissesthread, 5)
			self.runVomMisses(vonmissesthread, 7)
			self.runVomMisses(vonmissesthread, 10)

		self.runVomMissesGridFreq(2, ngridfreq, tgridfreq)
		self.runVomMissesGridFreq(3, ngridfreq, tgridfreq)
		self.runVomMissesGridFreq(5, ngridfreq, tgridfreq)
		self.runVomMissesGridFreq(7, ngridfreq, tgridfreq)
		self.runVomMissesGridFreq(10, ngridfreq, tgridfreq)


		self.scanVM(apfilename + ".vm2.resgf", 2)
		self.scanVM(apfilename + ".vm3.resgf", 3)
		self.scanVM(apfilename + ".vm5.resgf", 5)
		self.scanVM(apfilename + ".vm7.resgf", 7)
		self.scanVM(apfilename + ".vm10.resgf", 10)

	def fullAnalysisLoadAP2(self, ap2filename, analyzevm=-1, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800):
		self.loadnormalizedAP2(ap2filename)
		self.freqmin=float(freqmin)
		self.freqmax=float(freqmax)
		self.vmnumax=float(vmnumax)
		self.scanLS(self.freqmin, self.freqmax, 0, 14)

		if analyzevm == 1:
			self.runVomMisses(vonmissesthread, 2)
			self.runVomMisses(vonmissesthread, 3)
			self.runVomMisses(vonmissesthread, 5)
			self.runVomMisses(vonmissesthread, 7)
			self.runVomMisses(vonmissesthread, 10)

			self.runVomMissesGridFreq(2, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(3, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(5, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(7, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(10, ngridfreq, tgridfreq)


			self.scanVM(apfilename + ".vm2.resgf", 2)
			self.scanVM(apfilename + ".vm3.resgf", 3)
			self.scanVM(apfilename + ".vm5.resgf", 5)
			self.scanVM(apfilename + ".vm7.resgf", 7)
			self.scanVM(apfilename + ".vm10.resgf", 10)


	def fullAnalysis3(self, apfilename, ranal=2, gasvalue=0.00054, analyzevm=-1, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800):
		self.normalizeAP3(apfilename, ranal, gasvalue)
		self.freqmin=float(freqmin)
		self.freqmax=float(freqmax)
		self.vmnumax=float(vmnumax)
		
		self.scanLS(self.freqmin, self.freqmax, 0, 18)
		#self.scanLS(-1, -1, 0, 19)

		if analyzevm == 1:
			self.runVomMisses(vonmissesthread, 0)
			self.runVomMisses(vonmissesthread, 1)
			self.runVomMisses(vonmissesthread, 2)
			self.runVomMisses(vonmissesthread, 3)
			self.runVomMisses(vonmissesthread, 8)
			
			self.runVomMissesGridFreq(0, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(1, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(2, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(3, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(8, ngridfreq, tgridfreq)
			
			
			self.scanVM(apfilename + ".vm0.resgf", 0)
			self.scanVM(apfilename + ".vm1.resgf", 1)
			self.scanVM(apfilename + ".vm2.resgf", 2)
			self.scanVM(apfilename + ".vm3.resgf", 3)
			self.scanVM(apfilename + ".vm8.resgf", 8)

	def fullAnalysisLoadAP3(self, ap3filename, analyzevm=-1, vonmissesthread=48, freqmin=0.5e-06, freqmax=5.0e-06, vmnumax=100, ngridfreq=1000, tgridfreq=10800):
		self.loadnormalizedAP3(ap3filename)
		self.freqmin=float(freqmin)
		self.freqmax=float(freqmax)
		self.vmnumax=float(vmnumax)
		self.scanLS(self.freqmin, self.freqmax, 0, 18)
		
		if analyzevm == 1:
			self.runVomMisses(vonmissesthread, 0)
			self.runVomMisses(vonmissesthread, 1)
			self.runVomMisses(vonmissesthread, 2)
			self.runVomMisses(vonmissesthread, 3)
			self.runVomMisses(vonmissesthread, 8)
			
			self.runVomMissesGridFreq(0, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(1, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(2, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(3, ngridfreq, tgridfreq)
			self.runVomMissesGridFreq(8, ngridfreq, tgridfreq)
			
			
			self.scanVM(ap3filename + ".vm0.resgf", 0)
			self.scanVM(ap3filename + ".vm1.resgf", 1)
			self.scanVM(ap3filename + ".vm2.resgf", 2)
			self.scanVM(ap3filename + ".vm3.resgf", 3)
			self.scanVM(ap3filename + ".vm8.resgf", 8)
