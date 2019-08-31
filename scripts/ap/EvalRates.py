import string, os, sys
import numpy as np
import sys
import math
from scipy.stats import norm

#e = EvalRates()
#e.calculateRateWithoutExp(1, 0.7, 0e-08, 1.2e-4, 0.76, 1.75, 400)

class EvalRates:

	def __init__(self):
		return

	def getInstrumentPSF(self, instrumentID = 0, emin=100):
		#AGILE
		psf = 0
		if instrumentID == 0:
			if emin < 1700.: psf = 0.3 # deg.
			if emin < 1000.: psf = 1. # deg.
			if emin < 400.: psf = 2.4 # deg.calibration with Vela with ranal=2, 3, 4
			if emin < 100.: psf = 5.0 # deg. da provare

		return psf
	
	def getFluxScaleFactor(self, verbose=0, ranal= 2, emin = 100., emax = 10000., instrumentID = 0):
		ranal = float(ranal)
		emin = float(emin)
		emax = float(emax)
		if verbose == 1:
			print('ranal %.2f'%ranal)
			print('emin %d'% emin)
			print('emax %d'% emax)
		psf = self.getInstrumentPSF(instrumentID, emin)
		if verbose == 1:
			print('selected psf: ' + str(psf))
			print('selected ranal: ' + str(ranal))
		
		#to take into account the extension of the PSF
		fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranal))
		
		print('Fluxscalefactor based on PSF: %.4f' % fluxscalefactor)

		# ON - PSF region
		omega_ranal = 2.*np.pi*(1. - np.cos(ranal*(np.pi/180.))) #[sr]
		#AC = Andrew
		omega_ranalAC = (np.pi/180. * ranal)**2 * np.sin(np.pi/180. * 30) / (np.pi/180. * 30) #[sr]
	
		if verbose == 1:
			print('[sr]   - ' + str(omega_ranal))
			print('[sr] AC - ' + str(omega_ranalAC))
			
		return fluxscalefactor;

	def calculateRateWithoutExp(self, verbose = 0, ranal= -1, fluxsource = 0e-08, gasvalue=-1, gal = 0.7, iso = 10., emin = 100., emax = 50000., instrumentID = 0):
		fluxsource = float(fluxsource)
		gasvalue = float(gasvalue) #andare direttamente nei file .disp.conv.sky.gz del modello e prendere il valore da li'
		gal = float(gal)
		iso = float(iso)
		ranal = float(ranal)
		emin = float(emin)
		emax = float(emax)
		if verbose == 1:
			print('gasvalue %.5f'%gasvalue)
			print('gascoeff %.2f'% gal)
			print('isocoeff %.2f'% iso)
			print('emin %d'% emin)
			print('emax %d'% emax)
		psf = self.getInstrumentPSF(instrumentID, emin)
		if verbose == 1:
			print('selected psf: ' + str(psf))
			print('selected ranal: ' + str(ranal))

		#to take into account the extension of the PSF
		fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranal))

		if verbose == 1:
			print('fluxscalefactor based on PSF: ' + str(fluxscalefactor))

		# ON - PSF region
		omega_ranal = 2.*np.pi*(1. - np.cos(ranal*(np.pi/180.))) #[sr]
		#AC = Andrew
		omega_ranalAC = (np.pi/180. * ranal)**2 * np.sin(np.pi/180. * 30) / (np.pi/180. * 30) #[sr]

		if verbose == 1:
			print('[sr]   - ' + str(omega_ranal))
			print('[sr] AC - ' + str(omega_ranalAC))

		bkgdata = gasvalue*gal + iso*(10.**(-5)) # cts / [cm2 s sr]
		bkg_ON = bkgdata*omega_ranal # [cts] / [cm^2 s sr] * [sr] = [cts] / [cm^2 s]

		ctsgal = gasvalue*gal * omega_ranal #[cts] / [cm^2 s]
		ctsiso = iso*(10.**(-5)) * omega_ranal #[cts] / [cm^2 s]

		if verbose == 1:
			print('--------------')
			print('[cts] / [cm^2 s] - bkg_ON = ctsgal + ctsiso = ' + str(bkg_ON) + ' = ' + str(ctsgal) + ' + ' + str(ctsiso))

		#B) flux source constant
		fluxsource = fluxsource * fluxscalefactor #[cts / cm2 s]

		#Calculation of counts
		src_ON =  fluxsource #  [cts] / [cm^2 s]

		ctstot = bkg_ON + src_ON # [cts] / [cm^2 s]

		if verbose == 1:
			print('[cts] / [cm^2 s] - src_ON = ' + str(src_ON))
			print('[cts] / [cm^2 s] - ctstot (bkg + src) = ' + str(ctstot))
			print('-------------- Moltiplica il valore sopra per exp in [cm2 s]')

		return bkg_ON, src_ON

	#exposure [cm2 s]
	#fluxsource [cts] / [cm2 s]
	#gasvalue [cts] / [cm2 s sr]

	def calculateRateAndSNR(self, verbose = 0, ranal= -1, exposure = 40000, fluxsource = 0e-08, gasvalue=-1, gal = 0.7, iso = 10., emin = 100., emax = 50000., instrumentID = 0):
		bkg_ON, src_ON = self.calculateRateWithoutExp(verbose, ranal, fluxsource, gasvalue, gal, iso, emin, emax, instrumentID)
		ctstot = (bkg_ON + src_ON) * exposure # [cts]
		snr = src_ON * exposure / math.sqrt(float(ctstot))

		print('src_ON = %.3e'% src_ON)
		print('bkg_ON = %.3e'% bkg_ON)
		print('ctstot = %.3f'% ctstot)
		print('SNR = %.3f' % snr)
		print('--------------')

	def determinebestSNR(self, verbose=0, ranalstart= 0.1, ranalend =4.0, exposure = 40000, fluxsource = 0e-08, gasvalue=0.0006, gal = 0.7, iso = 10., emin = 100., emax = 10000., instrumentID = 0):

		for ranal in np.arange(ranalstart, ranalend, 0.1):
			bkg_ON, src_ON = self.calculateRateWithoutExp(verbose, ranal, fluxsource, gasvalue, gal, iso, emin, emax, instrumentID)
			ctstot = (bkg_ON + src_ON) * exposure # [cts]
			snr = src_ON * exposure / math.sqrt(float(ctstot))
			print('%.2f %.3f' % (ranal, snr))
		return
