import string, os, sys
import numpy as np
import sys
import math
from scipy.stats import norm
from astropy.io import fits
from matplotlib import gridspec
from EdpGrid import *

#e = EvalRates()
#e.calculateRateWithoutExp(1, 0.7, 0e-08, 1.2e-4, 0.76, 1.75, 400)

class MathUtils:
	def steradiansCone(this, ranal):
		return 2.*np.pi*(1. - np.cos(ranal*(np.pi/180.))) #[sr]
		
	def steradiansSquarePixel(this, dimpixel, theta=30):
		pixel = (np.pi/180. * dimpixel)**2
		return pixel * np.sin(np.pi/180. * theta) / (np.pi/180. * theta) #[sr]
		#omega_ranalAC = (np.pi/180. * ranal)**2 * np.sin(np.pi/180. * 30) / (np.pi/180. * 30) #[sr]
		
		
class PSFEval:

	def __init__(self):
		return

	# model functions
	def PowerLaw(self, x, slope):
		return x**(-slope)
	
	def IntPowerLaw(self, e1, e2, slope):
		return (((e2**(-slope + 1.))/(-slope + 1.)) - ((e1**(-slope + 1.))/(-slope + 1.)))

	
	def EvalPSFMean(self, emin=100, emax=10000, gindex=2.1, source_theta=30, source_phi=0):
		irf_path = os.environ['AGILE'] + '/model/scientific_analysis/data/'

		
		######## loading the IRF PSF and effective area
		psf_file = irf_path+"/AG_GRID_G0017_SFMG_H0025.psd.gz"
		aeff_file = irf_path+"/AG_GRID_G0017_SFMG_H0025.sar.gz"
		
		# reading the effective area
		hdulist_aeff = fits.open(aeff_file)
		
		wcs_tab_aeff = hdulist_aeff[2].data
		Energy_min = wcs_tab_aeff.field('ENERGY')
		theta = wcs_tab_aeff.field('POLAR_ANGLE')
		phi = wcs_tab_aeff.field('AZIMUTH_ANGLE')
		
		# computing the energy bin center
		Energy_min = Energy_min[0]
		Energy_max = [35, 50, 71, 100, 173, 300, 548, 1000, 1732, 3000, 5477, 10000, 20000, 50000]
		Energy = []
		for jene in range(len(Energy_min)):
			Energy.append(Energy_min[jene] + (Energy_max[jene] - Energy_min[jene])/2.)
		
		
		Energy = np.array(Energy)
		Energy_max = np.array(Energy_max)
		Energy_min = np.array(Energy_min)
		
		# Computing the PSF 
		theta_sel = source_theta
		phi_sel = source_phi
		
		theta = theta[0]
		phi = phi[0]
		
		hdulist_psf = fits.open(psf_file)
		
		wcs_tab_psf = hdulist_psf[2].data
		Rho = wcs_tab_psf.field('RHO')
		Rho = Rho[0]
		Psi = wcs_tab_psf.field('PSI')
		Psi = Psi[0]
		
		primary_psf = hdulist_psf[0].data
		PSF = np.zeros(len(Energy))
		for jphi in range(len(phi)):
			if (phi[jphi] == phi_sel):
				primary_psf_phi = primary_psf[jphi]
				#print "primary_psf_phi ", primary_psf_phi
				for jtheta in range(len(theta)):
					if (theta[jtheta] == theta_sel):
						primary_psf_phi_theta = primary_psf_phi[jtheta]
						for jene in range(len(Energy)):
							primary_psf_phi_theta_ene = primary_psf_phi_theta[jene]
							for jpsi in range(len(Psi)):
								if (Psi[jpsi] == 0):
									primary_psf_phi_theta_ene_psi = primary_psf_phi_theta_ene[jpsi]
									counts = []
									for jrho in range(len(Rho)):
										sph_annulus = (2.*math.pi*(math.cos((Rho[jrho]-0.05)*(math.pi/180.)) - math.cos((Rho[jrho]+0.05)*(math.pi/180.)))) #sr
										counts.append(primary_psf_phi_theta_ene_psi[jrho]*sph_annulus)
										
									total_counts = np.sum(counts)
									cr_counts = 0.68*total_counts
									radial_counts = 0
									for jrho in range(len(Rho)):
										radial_counts = radial_counts + counts[jrho]
										if (radial_counts >= cr_counts): 
											PSF[jene] = Rho[jrho]
											break
		
		
		####### PSF averaged over the source energy distribution
		
		
		# select energy range
		where_band = np.where((Energy_min >= emin) & (Energy_max <= emax))
		energy_band = Energy[where_band]
		PSF_band = PSF[where_band]
		energymin_band = Energy_min[where_band]
		energymax_band = Energy_max[where_band]
		
		int_source_band = self.IntPowerLaw(emin, emax, gindex)
		fract_source = []
		for je in range(len(energy_band)):
			int_source_bin = self.IntPowerLaw(energymin_band[je], energymax_band[je], gindex)
			fract_source.append(int_source_bin/int_source_band)
		
		PSF_norm = np.average(PSF_band, weights = fract_source)
		
		print('###########################################')
		print('#              SOURCE PSF                 #')
		print('###########################################')
		print('# - Energy min [MeV] = %d'% emin)
		print('# - Energy max [MeV] = %d'% emax)
		print('# - photon index [E^-gamma] = %.2f'% gindex)
		print('# - source off-axis [deg.] = %d'% theta_sel)
		print('# - source azimuthal angle [deg.] = %d'% phi_sel)
		print('# - normalized PSF [deg.] = %.4f'% PSF_norm)
		
		return PSF_norm
		


class EvalRates:

	def __init__(self):
		return

	def getInstrumentPSF(self, instrumentID = 0, gindex=2.1, emin=100, emax=10000, source_theta=30):
		#AGILE
		psf = 0
		psfc = PSFEval()
		if instrumentID == 0:
			if emin < 1700.: psf = 0.3 # deg.
			if emin < 1000.: psf = 1. # deg.
			if emin < 400.: psf = 2.3 # deg.
			if emin < 100.: psf = 5.0 # deg. da provare
		if instrumentID == 0:
			psf = psfc.EvalPSFMean(emin, emax, gindex, source_theta)

		return psf
	
	def getFluxScaleFactor(self, verbose=0,  gindex=2.1, ranal= 2, emin = 100., emax = 10000., instrumentID = 0, source_theta=30):
		ranal = float(ranal)
		emin = float(emin)
		emax = float(emax)
		if verbose == 1:
			print('ranal %.2f'%ranal)
			print('emin [MeV]: %d'% emin)
			print('emax [MeV]: %d'% emax)
		psf = self.getInstrumentPSF(instrumentID, gindex, emin, emax, source_theta)
		if verbose == 1:
			print('selected psf  : ' + str(psf))
			print('selected ranal: ' + str(ranal))
		
		#to take into account the extension of the PSF
		fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranal))
		print('Fluxscalefactor based on PSF : %.4f' % fluxscalefactor)
		
		#to take into account the spectral shape of the source (deviation from spectral index=2.1)
		edpGrid = EdpGrid()
		edp_file = os.environ['AGILE']+"/model/scientific_analysis/data/AG_GRID_G0017_SFMG_H0025.edp.gz"
		edpGrid.readData(edp_file)
		corrsi = edpGrid.detCorrectionSpectraFactorSimple(edpGrid, emin, emax, gindex)
		print('Correction spectra factor: %.2f'%corrsi)
		fluxscalefactor = fluxscalefactor / corrsi
		
		print('Fluxscalefactor based on PSF and correction spectra factor: %.4f' % fluxscalefactor)
		
		# ON - PSF region
		#omega_ranal = 2.*np.pi*(1. - np.cos(ranal*(np.pi/180.))) #[sr]
		mu = MathUtils()
		omega_ranal = mu.steradiansCone(ranal)
		
		if verbose == 1:
			print('omega    [sr]  : ' + str(omega_ranal))
			#print('omega AC [sr]  : ' + str(omega_ranalAC))
			
		return fluxscalefactor;

	def calculateRateWithoutExp(self, verbose = 0, ranal= -1, fluxsource = 0e-08, gasvalue=-1, gal = 0.7, iso = 10., emin = 100., emax = 10000., gindex=2.1, source_theta=30, instrumentID = 0):
		fluxsource = float(fluxsource)
		gasvalue = float(gasvalue) #andare direttamente nei file .disp.conv.sky.gz del modello e prendere il valore da li'
		gal = float(gal)
		iso = float(iso)
		ranal = float(ranal)
		emin = float(emin)
		emax = float(emax)
		if verbose == 1:
			print('gasvalue    %.5f'%gasvalue)
			print('gascoeff    %.2f'% gal)
			print('isocoeff    %.2f'% iso)
			print('emin [MeV]: %d'% emin)
			print('emax [MeV]: %d'% emax)
		psf = self.getInstrumentPSF(instrumentID, gindex, emin, emax, source_theta)
		if verbose == 1:
			print('selected psf   [deg]: ' + str(psf))
			print('selected ranal [deg]: ' + str(ranal))

		#to take into account the extension of the PSF
		fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranal))

		if verbose == 1:
			print('fluxscalefactor based on PSF: ' + str(fluxscalefactor))

		# ON - PSF region - steradians of a cone of angle \theta
		#omega_ranal = 2.*np.pi*(1. - np.cos(ranal*(np.pi/180.))) #[sr]
		mu = MathUtils()
		omega_ranal = mu.steradiansCone(ranal)

		if verbose == 1:
			print('[sr]    : ' + str(omega_ranal))
			#print('[sr] AC : ' + str(omega_ranalAC))

		bkgdata = gasvalue*gal + iso*(10.**(-5)) # cts / [cm2 s sr]
		print("absolute background (gasvalue * gal + iso*(10^-5)) [cts / cm2 s sr]: %3f"%bkgdata)
		bkg_ON = bkgdata*omega_ranal # [cts] / [cm^2 s sr] * [sr] = [cts] / [cm^2 s]

		ctsgal = gasvalue * gal * omega_ranal #[cts] / [cm^2 s]
		ctsiso = iso*(10.**(-5)) * omega_ranal #[cts] / [cm^2 s]

		if verbose == 1:
			print('--------------')
			print('bkg_ON = ctsgal + ctsiso [cts / cm^2 s]: ' + str(bkg_ON) + ' = ' + str(ctsgal) + ' + ' + str(ctsiso))

		#B) flux source constant
		fluxsource = fluxsource * fluxscalefactor #[cts / cm2 s]

		#Calculation of counts
		src_ON =  fluxsource #  [cts] / [cm^2 s]

		ctstot = bkg_ON + src_ON # [cts] / [cm^2 s]

		if verbose == 1:
			print('bkg_ON [cts / cm^2 s]:            %.3e' % bkg_ON)
			print('src_ON [cts / cm^2 s]:            %.3e' % src_ON)
			print('ctstot (bkg + src) [cts / cm2 s]: %.3e' % ctstot)
			#print('-------------- Moltiplica il valore sopra per exp in [cm2 s]')

		return bkg_ON, src_ON

	#exposure [cm2 s]
	#fluxsource [cts] / [cm2 s]
	#gasvalue [cts] / [cm2 s sr]

	def calculateRateAndSNR(self, verbose = 0, ranal= -1, exposure = 40000, fluxsource = 0e-08, gasvalue=-1, gal = 0.7, iso = 10., emin = 100., emax = 10000., gindex=2.1, source_theta=30, instrumentID = 0):
		bkg_ON, src_ON = self.calculateRateWithoutExp(verbose, ranal, fluxsource, gasvalue, gal, iso, emin, emax, gindex, source_theta, instrumentID)
		ctstot = (bkg_ON + src_ON) * exposure # [cts]
		snr = src_ON * exposure / math.sqrt(float(ctstot))

		print('src_ON [cts / cm^2 s]: %.3e'% src_ON)
		print('bkg_ON [cts / cm^2 s]: %.3e'% bkg_ON)
		print('ctstot [cts]:          %.3f'% ctstot)
		print('SNR: %.3f' % snr)
		print('--------------')
		
		alpha=1
		N_off = bkg_ON * exposure
		N_on  = src_ON * exposure
		A_part = ((1. + alpha)/alpha)*(N_on/(N_on + N_off))
		B_part = (1. + alpha)*(N_off/(N_on + N_off))
		S_lima = np.sqrt(2)*math.sqrt((N_on*np.log(A_part)) + (N_off*np.log(B_part)))
		print('S_lima: %.3f' % snr)

	def determinebestSNR(self, verbose=0, ranalstart= 0.1, ranalend =4.0, exposure = 40000, fluxsource = 0e-08, gasvalue=0.0006, gal = 0.7, iso = 10., emin = 100., emax = 10000., gindex=2.1, source_theta=30, instrumentID = 0):

		for ranal in np.arange(ranalstart, ranalend, 0.1):
			bkg_ON, src_ON = self.calculateRateWithoutExp(verbose, ranal, fluxsource, gasvalue, gal, iso, emin, emax, gindex, source_theta, instrumentID)
			ctstot = (bkg_ON + src_ON) * exposure # [cts]
			snr = src_ON * exposure / math.sqrt(float(ctstot))
			print('%.2f %.3f' % (ranal, snr))
		return
