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
		#AA
		sum1 = 0.0
		sum2 = 0.0
		sum3 = 0.0
		diml = len(expdataA)
		flux = np.zeros(diml)

		n=0
		for e_i in expdataA:
			if len(ctssimA) > 0:
				nb_i = ctssimA[n]
			else:
				nb_i = ctsdataA[n]
			n_i = ctsdataA[n]

			flux_i = nb_i / e_i
			#method 1
			s_ip = 0.5 + np.sqrt(nb_i + 0.25)
			s_im = -0.5 + np.sqrt(nb_i + 0.25)
			s_irms = np.sqrt(s_ip*s_ip + s_im*s_im) / 2.0
			s_i = s_irms / e_i

			sum1 += flux_i / (s_i*s_i)
			sum2 += 1.0 / (s_i*s_i)
			sum3 += flux_i

			flux_i = n_i / e_i
			flux[n] = flux_i

			n=n+1

		flux_mean = sum3 / n
		fluxw_mean = sum1 / sum2

		n=0
		for flux_i in flux:
			fluxpred_i = expdataA[n] * flux_mean
			s_i = np.sqrt(fluxpred_i) / expdataA[n]
			#aa_fluxw << (flux_i - fluxw_mean) / s_i
			if s_i != 0.0:
				fluxAA[n] = (flux_i - fluxw_mean) / (s_i*s_i)
			else:
				fluxAA[n] = 0.0

			n = n + 1

		return

	def normalizeAB(self, expdataA, ctsdataA, fluxAB1, fluxAB2, fluxAB3, ctssimA):
		sum1 = 0.0
		sum2 = 0.0
		sum3 = 0.0
		sum4 = 0.0
		sum5 = 0.0
		diml = len(expdataA)
		flux = np.zeros(diml)

		ab1_fluxw = []
		ab2_fluxw = []
		ab3_fluxw = []
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
			sum3 += flux_i
			sum4 += nb_i
			sum5 += e_i

			flux_i = n_i / e_i
			flux[n] = flux_i
			n = n + 1

		fluxw_mean1 = sum1 / sum2
		fluxw_mean2 = sum3 / n
		fluxw_mean3 = sum4 / sum5
		n=0
		for flux_i in flux:
			e_i = expdataA[n]

			fluxpred_i = e_i * fluxw_mean1
			s_i = np.sqrt(fluxpred_i) / e_i
			#ab1_fluxw << (flux_i - fluxw_mean1) / s_i
			fluxAB1[n] = (flux_i - fluxw_mean1) / (s_i*s_i)

			fluxpred_i = e_i * fluxw_mean2
			s_i = np.sqrt(fluxpred_i) / e_i
			#ab2_fluxw << (flux_i - fluxw_mean2) / s_i
			fluxAB2[n] = (flux_i - fluxw_mean2) / (s_i*s_i)

			#red text method
			fluxpred_i = e_i * fluxw_mean3
			s_i = np.sqrt(fluxpred_i) / e_i
			#ab3_fluxw << (flux_i - fluxw_mean3) / s_i
			fluxAB3[n] = (flux_i - fluxw_mean3) / (s_i*s_i)

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
		#0 flux cts/exp
		#1 variance of the flux
		#2 def_fluxw
		#3 aa_fluxw
		#4 ab1_fluxw
		#5 ab2_fluxw
		#6 ab3_fluxw
		#7 def_fluxwB
		#8 aa_fluxwB
		#9 ab1_fluxwB
		#10 ab2_fluxwB
		#11 ab3_fluxwB
		#12 cts
		#13 exp
		self.res = []
		return

	def loadDataAPAGILE(self, apfile):
		print('Loading data...')
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

	def calculateCtsFromConstantModel(self, ranal= 1, gasvalue=0, gal = 0.7, iso = 20., emin = 100., emax = 10000., instrumentID = 0):
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

	def normalizeAP(self, apfile):

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

		#0 flux cts/exp
		#1 variance of the flux
		#2 def_fluxw
		nap.normalizeDef(self.expdataA, self.ctsdataA, self.res[:,2], self.res[:,0], self.res[:,1], [])

		#3 aa_fluxw
		nap.normalizeAA(self.expdataA, self.ctsdataA, self.res[:,3], [])

		#4 ab1_fluxw
		#5 ab2_fluxw
		#6 ab3_fluxw
		nap.normalizeAB(self.expdataA, self.ctsdataA, self.res[:,4], self.res[:,5], self.res[:,6], [])

		if len(self.ctssimA) == 0:
			self.calculateCtsFromConstantModel(2, 0.0006, 1, 30, 100, 10000, 0)

		if len(self.ctssimA) > 0:
			#0 flux cts/exp
			#1 variance of the flux
			#7 def_fluxwB
			nap.normalizeDef(self.expdataA, self.ctsdataA, self.res[:,7], self.res[:,0], self.res[:,1], self.ctssimA)
			#8 aa_fluxwB
			nap.normalizeAA(self.expdataA, self.ctsdataA, self.res[:,8], self.ctssimA)
			#9 ab1_fluxwB
			#10 ab2_fluxwB
			#11 ab3_fluxwB
			nap.normalizeAB(self.expdataA, self.ctsdataA, self.res[:,9], self.res[:,10], self.res[:,11], self.ctssimA)

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
		print("* Write ap2 file: tstart tstop exp[cm2 s] cts 0:flux[cts/exp] 1:flux_var 2:def_fluxw 3:aa_fluxw 4:ab1_fluxw 5:ab2_fluxw 6:ab3_fluxw 7:def_fluxwB 8:aa_fluxwB 9:ab1_fluxwB 10:ab2_fluxwB 11:ab3_fluxwB")
		print("  B means normalization for a background model: calculateCtsFromConstantModel(2, 0.0006, 1, 30, 100, 10000, 0)")

		self.writeVonMisses(apfile, 2)
		self.writeVonMisses(apfile, 3)
		self.writeVonMisses(apfile, 4)
		self.writeVonMisses(apfile, 10)
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
		print("* Write Von Misses file: tcenter  "+str(ii)+":flux 1:flux_var")

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

		return pls, power.max(), maxf

	def scanLS(self, minfreq=-1, maxfreq=-1):
		fileclean = open(self.apfile + ".apres","w")
		for i in range(0,14):
			pls, pmax, maxf = self.calculateLS(0 ,0, i, minfreq, maxfreq)
			daymax = 1 / maxf / 86400.
			print("LS " + str(i) + ' ' + str(pls) + ' ' + str(pmax) + ' ' + str(maxf) + ' ' + str(daymax))
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
						#print(val)
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

		if plot == 1:
			plt.plot(frequency, power)
			plt.gca().set_xscale("log")
			#plt.gca().set_yscale("log")
			plt.show()

		return maxf, maxp

	def scanVM(self, vmfilename):
		vmtype = vmfilename.split(".ap.")[1].split(".")[0]
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
		cmd = "module load icc-18.0.1; "+os.environ['AGILE']+"/bin/eval_vonmises.prg "+str(nthreads)+" 0 0.5e-06 5.0e-06 0 100 < " + apfilename + ".vm"+str(ii)+" > " + apfilename + ".vm"+str(ii)+".res"
		os.system(cmd)
		cmd = "module load icc-18.0.1; "+os.environ['AGILE']+"/bin/grid_freq.prg 0.5e-06 5.0e-06 1000 600 < " + apfilename + ".vm"+str(ii)+".res > " + apfilename + ".vm"+str(ii)+".resgf"
		os.system(cmd)

	def fullAnalysis(self, apfilename, analyzevm=-1, vonmissesthread=48):
		self.normalizeAP(apfilename)
		if analyzevm == 1:
			self.runVomMisses(vonmissesthread, 2)
			self.runVomMisses(vonmissesthread, 3)
			self.runVomMisses(vonmissesthread, 5)
			#runVomMisses(48, 10)

		self.scanLS(0.5e-6, 5e-6)
		self.scanVM(apfilename + ".vm2.resgf")
		self.scanVM(apfilename + ".vm3.resgf")
		self.scanVM(apfilename + ".vm5.resgf")
		#self.scanVM(apfilename + ".vm10.resgf")
