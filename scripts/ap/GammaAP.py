import scipy.signal as signal
from astropy.stats import LombScargle
import string, os, sys
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from EvalRates import *

class NormalizeAP:

	def __init__(self):
		return

	# https://docs.google.com/document/d/1FLSHiYUU1Lcl2JFkwAMKMU0GG9GbX4Ou0jUADRYIPr8/edit?usp=sharing
	# Section 2
	#
	# [1] Robin Corbet, LAT Light Curve Analysis: Aperture Photometry and Periodicity Searches, slides at COSPAR CBW 2010
	# [2] Gehrels, 1986, ApJ, 303, 336
	# [3] Kraft, Burrows, & Nousek, 1991, ApJ, 374, 344

	#def_fluxw
	# Variance weighting corrected for low count statistic - Weighted Power Spectrum	 [1] [2]
	def normalizeDef(self, expdataA, ctsdataA, fluxw, flux, s, ctssimA):
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

			flux_i = nb_i / e_i
			s_ip = 0.5 + np.sqrt(nb_i + 0.25)
			s_im = -0.5 + np.sqrt(nb_i + 0.25)
			s_irms = np.sqrt(s_ip*s_ip + s_im*s_im) / 2.0
			s_i = s_irms / e_i
			sum1 += flux_i / (s_i*s_i)
			sum2 += 1.0 / (s_i*s_i)

			flux_i = n_i / e_i
			flux[n] = flux_i
			s[n] = s_i

			n = n + 1

		fluxw_mean = sum1 / sum2

		n=0
		for flux_i in flux:
			s_i =  s[n]
			if s_i != 0.0:
				#def_fluxw << (flux_i - fluxw_mean) / s_i
				fluxw[n] = (flux_i - fluxw_mean) / (s_i*s_i)
			else:
				fluxw[n] = 0.0

			n = n + 1

		return

	def normalizeAA(self, expdataA, ctsdataA, fluxAA, ctssimA):
		return

	def normalizeAB(self, expdataA, ctsdataA, fluxAB1, fluxAB2, fluxAB3, ctssimA):
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
		#0 flux cts/exp
		#1 variance of the flux
		#2 def_fluxw
		#3 aa_fluxw
		#4 ab1_fluxw
		#5 ab2_fluxw
		#6 ab3_fluxw
		#7 def_fluxwS
		#8 aa_fluxwS
		#9 ab1_fluxwS
		#10 ab2_fluxwS
		#11 ab3_fluxwS
		self.res = []
		return

	def loadDataAPAGILE(self, apfile):
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

	def calculateCtsFromConstantModel(ranal= 1, fluxsource = 0e-08, gasvalue=0, gal = 0.7, iso = 20., emin = 100., emax = 10000., instrumentID = 0):
		if self.diml == 0:
			print('Load the data')
			return

		rate = EvalRates()
		ctstotrate = rate.calculateRateWithoutExp(0, ranal, 0e-08,  gasvalue, gal, iso, emin, emax, instrumentID)
		self.ctssimA = np.zeros(self.diml)
		ii = 0
		for exp in self.expdataA:
			ctsdatabkg = np.random.poisson(ctstotrate * exp)
			self.ctssimA[ii]
			ii += 1
		return

	def normalizeAP(self, apfile):

		if self.diml == 0:
			print('Load the data')
			self.loadDataAPAGILE(apfile)
			#return


		self.res = np.zeros((self.diml, 12))
		nap = NormalizeAP()

		#0 flux cts/exp
		#1 variance of the flux
		#2 def_fluxw
		nap.normalizeDef(self.expdataA, self.ctsdataA, self.res[:,2], self.res[:,0], self.res[:,1], [])

		#3 aa_fluxw
	#	nap.normalizeAA(self.expdataA, self.ctsdataA, self.res[:,3], [])

		#4 ab1_fluxw
		#5 ab2_fluxw
		#6 ab3_fluxw
	#	nap.normalizeAB(self.expdataA, self.ctsdataA, self.res[:,4], self.res[:,5], self.res[:,6], [])

		if len(self.ctssimA) > 0:
			#0 flux cts/exp
			#1 variance of the flux
			#7 def_fluxwS
			nap.normalizeDef(self.expdataA, self.ctsdataA, self.res[:,7], self.res[:,0], self.res[:,1], self.ctssimA)
			#8 aa_fluxwS
	#		nap.normalizAA(self.expdataA, self.ctsdataA, self.res[:,8], self.ctssimA)
			#9 ab1_fluxwS
			#10 ab2_fluxwS
			#11 ab3_fluxwS
	#		nap.normalizeAB(self.expdataA, self.ctsdataA, self.res[:,9], self.res[:,10], self.res[:,11], self.ctssimA)

		fileclean = open(apfile + ".apB","w")

		n = 0
		for x in self.expdataA:
			fileclean.write(str(self.tstartA[n]) + " " + str(self.tstopA[n]) + " " + str(self.expdataA[n]) + " " + str(int(self.ctsdataA[n])) + " " + str(self.res[n,0]) + " " + str(self.res[n,1]) + " " + str(self.res[n,2]) + "\n")
			n = n + 1

		fileclean.close()
		print('Load the data2')
		return

	def calculateLS(self, verbose=0, plot=1):
		#normalization standard model log psd
		ls = LombScargle(self.tstartA, self.ctsdataA, normalization='standard')
		frequency, power = ls.autopower(minimum_frequency=1e-7, maximum_frequency=10e-3, samples_per_peak=100)
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

		return pls, maxf
