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
import sys
import math
from scipy.stats import norm
from astropy.io import fits
from matplotlib import gridspec
from Edp import *
from PSFEval import *
from APSignificance import *
from MathUtils import *

#e = EvalRates()
#e.calculateRateWithoutExp(1, 0.7, 0e-08, 1.2e-4, 0.76, 1.75, 400)

class EvalRates:

	def __init__(self):
		return

	##########################################################################
	def getInstrumentPSF(self, instrumentID = 0, gindex=2.1, emin=100, emax=10000, source_theta=30, verbose=0):
		#AGILE
		psf = 0
		psfc = PSFEval()
#		if instrumentID == 0:
#			if emin < 1700.: psf = 0.3 # deg.
#			if emin < 1000.: psf = 1. # deg.
#			if emin < 400.: psf = 2.3 # deg.
#			if emin < 100.: psf = 5.0 # deg. da provare
		if instrumentID == 0:
			psf = psfc.EvalPSFMean(emin=emin, emax=emax, gindex=gindex, source_theta=source_theta, verbose=verbose)

		return psf
	
	##########################################################################
	def getFluxScaleFactor(self, verbose=0,  gindex=2.1, ranal= 2, emin = 100., emax = 10000., instrumentID = 0, source_theta=30):
	
		if instrumentID > 0:
			return 0
	
		ranal = float(ranal)
		emin = float(emin)
		emax = float(emax)
		if verbose == 1:
			print('###########################################')
			print('#           getFluxScaleFactor            #')
			print('###########################################')
			print('input ranal %.2f'%ranal)
			print('input emin [MeV]: %d'% emin)
			print('input emax [MeV]: %d'% emax)
			print('input selected ranal: ' + str(ranal))
		
		#Integral PSF evaluation

		psf = self.getInstrumentPSF(instrumentID=instrumentID, gindex=gindex, emin=emin, emax=emax, source_theta=source_theta, verbose=verbose)
		if verbose == 1:
			print('selected psf  : %.4f' % psf)
		
		#fluxscalefactor to take into account the extension of the PSF and the fraction enclosed into the ranal
		#fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranal))
		#if verbose == 1:
		#	print('Fluxscalefactor based on PSF enclosed fraction: %.4f' % fluxscalefactor)
		
		#fluxscalefactor2
		psfc = PSFEval()
		fluxscalefactor = psfc.EvalPSFScaleFactor(ranal=ranal, emin=emin, emax=emax, gindex=gindex, source_theta=source_theta, verbose=verbose)
		if verbose == 1:
			print('Fluxscalefactor2 based on PSF enclosed fraction: %.4f' % fluxscalefactor)
		
		#fluxscalefactor to take into account the spectral shape of the source (deviation from spectral index=2.1 of calculated exposure)
		edp = Edp()
		corrsi = edp.detCorrectionSpectraFactorSimple(emin, emax, gindex)
		if verbose == 1:
			print('Fluxscalefactor based on expfluxcorrection (correction for expsoure spectra factor): %.2f'%corrsi)
		
		fluxscalefactor = fluxscalefactor / corrsi
		if verbose == 1:
			print('Total fluxscalefactor: %.4f' % fluxscalefactor)
		
		if verbose == 1:
			print('###########################################')
		
		# ON - PSF region
		#omega_ranal = 2.*np.pi*(1. - np.cos(ranal*(np.pi/180.))) #[sr]
		#mu = MathUtils()
		#omega_ranal = mu.steradiansCone(ranal)
		
		#if verbose == 1:
		#	print('omega    [sr]  : ' + str(omega_ranal))
		#	print('omega AC [sr]  : ' + str(omega_ranalAC))
			
		return fluxscalefactor;

	##########################################################################
	def calculateRateWithoutExp(self, verbose = 0, ranalS= -1, fluxsource = 0e-08, gasvalue=-1, gal = 0.7, iso = 10., emin = 100., emax = 10000., gindex=2.1, source_theta=30, instrumentID = 0):

		gasvalue = float(gasvalue) #andare direttamente nei file .disp.conv.sky.gz del modello e prendere il valore da li'

		if verbose == 1:
			print('#######################################################################')
			print('#                     calculateRateWithoutExp                         #')
			print('#######################################################################')
			print('input gasvalue    %.5f'%gasvalue)
			print('input gascoeff    %.2f'% gal)
			print('input isocoeff    %.2f'% iso)
			print('input emin [MeV]  %d'% emin)
			print('input emax [MeV]  %d'% emax)
			print('input ranalS      %d'% ranalS)
		#psf = self.getInstrumentPSF(instrumentID, gindex, emin, emax, source_theta)
		#if verbose == 1:
		#	print('selected psf   [deg]: ' + str(psf))
		#	print('selected ranalS [deg]: ' + str(ranalS))

		#to take into account the extension of the PSF
		#fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranalS))
		fluxscalefactor = self.getFluxScaleFactor(verbose = verbose, gindex=gindex, ranal= ranalS, emin = emin, emax = emax, instrumentID = instrumentID, source_theta=source_theta)
		
		if verbose == 1:
			print('S fluxscalefactor (PSF + expfluxcorrection):  %.3f' %fluxscalefactor)

		# ON - PSF region - steradians of a cone of angle \theta
		#omega_ranal = 2.*np.pi*(1. - np.cos(ranalS*(np.pi/180.))) #[sr]
		mu = MathUtils()
		omega_ranal = mu.steradiansCone(ranalS)

		if verbose == 1:
			print('[sr]    %.3f' %omega_ranal)
			#print('[sr] AC : ' + str(omega_ranalAC))

		rate_gal = gasvalue * gal * omega_ranal #[cts] / [cm^2 s]
		rate_iso = iso*(10.**(-5)) * omega_ranal #[cts] / [cm^2 s]

		bkgdata = gasvalue * gal + iso * (10.**(-5)) # cts / [cm2 s sr]
		#bkgdata = rate_gal + rate_iso
		if verbose == 1:
			print('B absolute background (gasvalue * gal + iso*(10^-5)) [cts / cm2 s sr]: %.3e'%bkgdata)
		
		rate_bkg_ON = bkgdata * omega_ranal # [cts] / [cm^2 s sr] * [sr] = [cts] / [cm^2 s]
		
		if verbose == 1:
			print('B rate_bkg_ON = rate_gal + rate_iso [cts / cm^2 s]: %.3e =  %.3e + %.3e' %( rate_bkg_ON, rate_gal, rate_iso))

		#Calculation of counts
		rate_src_ON =  fluxsource * fluxscalefactor #[cts / cm2 s]

		rate_tot = rate_bkg_ON + rate_src_ON # [cts] / [cm^2 s]
		
		flux_src_ON = fluxsource
		
		flux_tot = rate_bkg_ON + flux_src_ON
		
		if verbose == 1:
			print('rate_bkg_ON                          [cts / cm^2 s]:     %.3e' % rate_bkg_ON)
			print('rate_src_ON                          [cts / cm^2 s]:     %.3e' % rate_src_ON)
			print('rate_tot (bkg + src)                 [cts / cm^2 s]:     %.3e' % rate_tot)
			print('flux_src_ON                          [cts / cm^2 s]:     %.3e' % flux_src_ON)
			print('flux_tot (rate_bkg_ON + flux_src_ON) [cts / cm^2 s]:     %.3e' % flux_tot)
			#print('-------------- Moltiplica il valore sopra per exp in [cm2 s]')

		if verbose == 1:
			print('#######################################################################')
				
		return rate_bkg_ON, rate_src_ON

	##########################################################################
	#exposure = exposure in [cm2 s]
	#fluxsource = flux of S in [cts] / [cm2 s]
	#gasvalue = galactic diffuse emission coefficient (from SKY002) in [cts] / [cm2 s sr]
	#gal = gal coefficient
	#iso = iso coefficient
	#ranalS = radius of analysis of S
	#ranalB = radius of analysis of B. If the background is evaluated with AGILE MLE, ranalB=10, that is the usual radius of analysis used for the evalutation of gal and iso coefficients
	#alpha: alpha coefficient of Li&Ma (17) equation. If -1, let tool determine alpha
	def calculateRateAndSig(self, verbose = 0, ranalS= -1, exposure = 40000, fluxsource = 0e-08, gasvalue=-1, gal = 0.7, iso = 10., emin = 100., emax = 10000., gindex=2.1, source_theta=30, instrumentID = 0, ranalB = 10, alpha=-1):

		aps = APSignificance()
		
		rate_bkg_ON, rate_src_ON = self.calculateRateWithoutExp(verbose=verbose, ranalS=ranalS, fluxsource=fluxsource, gasvalue=gasvalue, gal=gal, iso=iso, emin=emin, emax=emax, gindex=gindex, source_theta=source_theta, instrumentID=instrumentID)
		
		#ctstot = (rate_bkg_ON + rate_src_ON) * exposure # [cts]
		
		ctsS = rate_src_ON * exposure
		ctsB = rate_bkg_ON * exposure
		
		snr = aps.SNR(verbose=verbose, ctsS=ctsS, ctsB=ctsB)
		
		ctstot = ctsS + ctsB
		N_off = rate_bkg_ON * exposure
		N_on  = (rate_bkg_ON + rate_src_ON) * exposure
		
		Slima = aps.lima(verbose=verbose, N_on = N_on, N_off = N_off, ranalS=ranalS)
		
		Sa = aps.Sa(verbose=verbose, ctsTOT=N_on, ctsB=N_off)
		
		if verbose >= 2:
			print('###########################################')
			print('#         calculateRate Sig UL Sens       #')
			print('###########################################')
			print('rate_bkg_ON [cts / cm^2 s]: %.3e'% rate_bkg_ON)
			print('rate_src_ON [cts / cm^2 s]: %.3e'% rate_src_ON)
			print('cts_tot_ON [cts]:           %.3f'% ctstot)
			print('ranalS:                     %.2f'% ranalS)
			print('N_ON  [cts]:                %.3f'% N_on)
			print('N_OFF [cts]:                %.3f'% N_off)
			print('Sig_SNR:                    %.3f'% snr)
			print('Sig_lima:                   %.3f'% Slima)
			print('Sig_a:                      %.3f'% Sa)
			
		aps = APSignificance()
		#upper limit
		fluxscalefactor = self.getFluxScaleFactor(verbose = verbose, gindex=gindex, ranal= ranalS, emin = emin, emax = emax, instrumentID = instrumentID, source_theta=source_theta)
		
		N_sourceUL, SignUL = self.calcCountsLimit(2, ctsB, ranalS, alpha=alpha)
		rateUL = (N_sourceUL / exposure)
		fluxUL = rateUL / fluxscalefactor
		
		if verbose >= 2:
			print('UL:                     ')
			print('cts_UL:                     %.1f'% N_sourceUL)
			print('Sign_UL:                    %.2f'% SignUL)
			print('rate_UL:                    %.2e'% rateUL)
			print('flux_UL:                    %.2e'% fluxUL)

		
		#sensitivity
		N_sourceSens, SignSens = self.calcCountsLimit(4, ctsB, ranalS, alpha=alpha)
		rateSens = (N_sourceSens / exposure)
		fluxSens = rateSens / fluxscalefactor
		#N_sourceUL, SignUL, rateUL, fluxUL %.1f %.2f %.2e %.2e
			
		if verbose >= 2:
			print('Sensitivity:                  ')
			print('cts_Sens:                   %.1f'% N_sourceSens)
			print('Sign_Sens:                  %.2f'% SignSens)
			print('rate_Sens:                  %.2e'% rateSens)
			print('flux_Sens:                  %.2e'% fluxSens)
		
		if verbose >= 1:
			print('###########################################')
		
		return	

	##########################################################################
	#ranalstartS = radius of analysis of S
	#ranalendS = radius of analysis of S
	#ranalB = radius of analysis of B. If the background is evaluated with AGILE MLE, ranalB=10, that is the usual radius of analysis used for the evalutation of gal and iso coefficients
	#alpha: alpha coefficient of Li&Ma (17) equation. If -1, let tool determine alpha
	def determinebestSig(self, verbose=0, ranalstartS= 0.1, ranalendS =4.0, exposure = 40000, fluxsource = 0e-08, gasvalue=0.0006, gal = 0.7, iso = 10., emin = 100., emax = 10000., gindex=2.1, source_theta=30, instrumentID = 0, ranalB = 10, alpha=-1):

		aps = APSignificance()
		
		for ranalS in np.arange(ranalstartS, ranalendS+0.1, 0.1):
			
			rate_bkg_ON, rate_src_ON = self.calculateRateWithoutExp(verbose=verbose, ranalS=ranalS, fluxsource=fluxsource, gasvalue=gasvalue, gal=gal, iso=iso, emin=emin, emax=emax, gindex=gindex, source_theta=source_theta, instrumentID=instrumentID)
			
			ctsS = rate_src_ON * exposure
			ctsB = rate_bkg_ON * exposure
			
			snr = aps.SNR(verbose=verbose, ctsS=ctsS, ctsB=ctsB)
			
			N_off = ctsB
			N_on  = ctsB + ctsS
			
			lima = aps.lima(verbose=verbose, N_on = N_on, N_off = N_off, ranalS=ranalS, alpha=alpha)
			
			Sa = aps.Sa(verbose=verbose, ctsTOT=N_on, ctsB=N_off)
			
			#upper limit
			fluxscalefactor = self.getFluxScaleFactor(verbose = verbose, gindex=gindex, ranal= ranalS, emin = emin, emax = emax, instrumentID = instrumentID, source_theta=source_theta)
			
			N_sourceUL, SignUL = self.calcCountsLimit(2, ctsB, ranalS, alpha=alpha)
			rateUL = (N_sourceUL / exposure)
			fluxUL = rateUL / fluxscalefactor
			
			#sensitivity
			N_sourceSens, SignSens = self.calcCountsLimit(4, ctsB, ranalS, alpha=alpha)
			rateSens = (N_sourceSens / exposure)
			fluxSens = rateSens / fluxscalefactor
			
			fluxSource = rate_src_ON / fluxscalefactor
			
			#saa = -1
			print('%.2f %.2f %.2f %.2f - %3d %.2e %.2e %3d %.2e - %.1f %.2f %.2e %.2e - %.1f %.2f %.2e %.2e' % (ranalS, snr, lima, Sa, ctsS, rate_src_ON, fluxSource, ctsB, rate_bkg_ON, N_sourceUL, SignUL, rateUL, fluxUL, N_sourceSens, SignSens, rateSens, fluxSens))
		return

	#ctsB
	def calcCountsLimit(self, Sign, ctsB, ranalS, alpha=-1, algorithm=2):
		
		aps = APSignificance()
		N_source = 0.1
		Ss = -1
		N_on = N_source + ctsB
		while(1):
			if (ctsB > 0.):
				N_on = N_source + ctsB
				Ss = -1
				if algorithm == 1:
					Ss = aps.lima(verbose=0, alpha=alpha, N_on = N_on, N_off = ctsB, ranalS=ranalS)
				if algorithm == 2:
					Ss = aps.Sa(verbose=0, ctsTOT=N_on, ctsB=ctsB)
				#print(N_source)
				if (Ss >= Sign):
					break
				else:
					N_source += 0.1
			else:
				N_source = 0.
				break

		return N_source, Ss
		#F_lim = (N_source/expo_on)/PSF_norm

		#print ("# Sensitivity parameters:")
		#print ("# - sigma = %.2f"% Sign[js])
		#print ("# - N_min = %.2f"% N_min[jmin])
		#print ("# Sensitivity [phot/cm2/s] = %.2e"% F_lim)
		#print ("###########################################")
					
		return

