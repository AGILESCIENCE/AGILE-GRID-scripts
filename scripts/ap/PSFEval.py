# DESCRIPTION
#       Agileap: AGILE Observatory Aperture Photometry Analysis
#       Main author of this file: Valentina Fioretti
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
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import norm
import matplotlib.pyplot as plt

class PSFEval:

	def __init__(self):
		return

	# model functions
	def PowerLaw(self, x, slope):
		return x**(-slope)

	def IntPowerLaw(self, e1, e2, slope):
		return (((e2**(-slope + 1.))/(-slope + 1.)) - ((e1**(-slope + 1.))/(-slope + 1.)))
	
	# King from FERMI
	def king_profile(self, x,sigma, gamma, B):
		return ((1./(2.*np.pi*(sigma**2)))*(1. - (1./gamma))*((1. + ((x**2)/(2.*(sigma**2)*gamma)))**(-gamma)))*B


	def EvalPSFMean(self, emin=100, emax=10000, gindex=2.1, source_theta=30, source_phi=0, verbose=0):
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
		
		if verbose == 1:
			print('###########################################')
			print('#              SOURCE PSF                 #')
			print('###########################################')
			print('# - Energy min [MeV] = %d'% emin)
			print('# - Energy max [MeV] = %d'% emax)
			print('# - photon index [E^-gamma] = %.2f'% gindex)
			print('# - source off-axis [deg.] = %d'% theta_sel)
			print('# - source azimuthal angle [deg.] = %d'% phi_sel)
			print('# - normalized PSF [deg.] = %.4f'% PSF_norm)
			print('###########################################')
		
		return PSF_norm


	def EvalPSFScaleFactor(self, ranal=4, emin=100, emax=10000, gindex=2.1, source_theta=30, source_phi=0, verbose=0):
	
		PSF_norm = self.EvalPSFMean(emin=emin, emax=emax, gindex=gindex, source_theta=source_theta, source_phi=source_phi, verbose=verbose)
		
		fluxscalefactor = math.fabs(1-2*norm(0,  PSF_norm).cdf(ranal))
		
		if verbose == 1:
			print('# - PSF scale factor [deg.] = %.4f'% fluxscalefactor)
		
		return fluxscalefactor
	
	
	def EvalPSFScaleFactor2(self, ranal=4, emin=100, emax=10000, gindex=2.1, source_theta=30, source_phi=0, verbose=0, plot_flag=False):
	
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

		fig = plt.figure(1,figsize=[7,6])
		ax = fig.add_subplot(111)
		
		if plot_flag:
			
			fig.tight_layout()
			ax.set_title(r'AGILE PSF ('+str(int(emin))+' < E < '+str(int(emax))+r', $\Gamma=$'+str(gindex)+')')


		####### PSF averaged over the source energy distribution

		# select energy range
		where_band = np.where((Energy_min >= emin) & (Energy_max <= emax))
		energy_band = Energy[where_band]
		#PSF_band = PSF[where_band]
		energymin_band = Energy_min[where_band]
		energymax_band = Energy_max[where_band]

		int_source_band = self.IntPowerLaw(emin, emax, gindex)
		fract_source = []
		for je in range(len(energy_band)):
			int_source_bin = self.IntPowerLaw(energymin_band[je], energymax_band[je], gindex)
			fract_source.append(int_source_bin/int_source_band)


		IRF_bin_rho = 0.1
		counter_eband = 0
		primary_psf = hdulist_psf[0].data
		density_PSF = np.zeros(len(Rho))
		counts_PSF = np.zeros(len(Rho))
		radius_PSF = np.zeros(len(Rho))
		for jphi in range(len(phi)):
			if (phi[jphi] == phi_sel):
				primary_psf_phi = primary_psf[jphi]
				for jtheta in range(len(theta)):
					if (theta[jtheta] == theta_sel):
						primary_psf_phi_theta = primary_psf_phi[jtheta]
						for jene in range(len(Energy)):
							primary_psf_phi_theta_ene = primary_psf_phi_theta[jene]
							for jpsi in range(len(Psi)):
								if (Psi[jpsi] == 0):
									primary_psf_phi_theta_ene_psi = primary_psf_phi_theta_ene[jpsi]
									if ((Energy[jene] >= emin) & (Energy[jene] <= emax)):
										modo_density_PSF = []
										for jrho in range(len(Rho)):
											radius_PSF[jrho] = Rho[jrho]
											density_PSF[jrho] = density_PSF[jrho] + primary_psf_phi_theta_ene_psi[jrho]*fract_source[counter_eband]
											modo_density_PSF.append(primary_psf_phi_theta_ene_psi[jrho])
											sph_annulus = (2.*np.pi*(np.cos((Rho[jrho]-IRF_bin_rho/2.)*(np.pi/180.)) - np.cos((Rho[jrho]+IRF_bin_rho/2.)*(np.pi/180.)))) #sr
											counts_PSF[jrho] = counts_PSF[jrho] + primary_psf_phi_theta_ene_psi[jrho]*fract_source[counter_eband]*sph_annulus
										counter_eband+=1
										# plotting the PSF for each energy bin
										modo_density_PSF = np.array(modo_density_PSF)
										top_mono = np.max(modo_density_PSF)
										ax.plot(radius_PSF, modo_density_PSF/top_mono, linewidth=1.5, label='E = '+str(Energy[jene])+' MeV')
									"""
									counts = []
									for jrho in range(len(Rho)):
										sph_annulus = (2.*np.pi*(np.cos((Rho[jrho]-IRF_bin_rho/2.)*(np.pi/180.)) - np.cos((Rho[jrho]+IRF_bin_rho/2.)*(np.pi/180.)))) #sr
										counts.append(primary_psf_phi_theta_ene_psi[jrho]*sph_annulus)
										
									total_counts = np.sum(counts)
									cr_counts = 0.68*total_counts
									radial_counts = 0
									for jrho in range(len(Rho)):
										radial_counts = radial_counts + counts[jrho]
										if (radial_counts >= cr_counts): 
											PSF[jene] = Rho[jrho]
											break
									"""

		# Fitting with a king profile

		# normalization
		density_PSF_norm = np.zeros(len(density_PSF))
		tot_rate = np.max(density_PSF)
		for jann in range(len(density_PSF_norm)):
			density_PSF_norm[jann] = density_PSF[jann]/tot_rate

				
		# king fit
		p, cov = curve_fit(self.king_profile, radius_PSF, density_PSF_norm, maxfev=1000000*(len(radius_PSF)+1))
					
		print('King fit result:')
		print(p)
		perr = np.sqrt(np.diag(cov))
		print('King fit result 1 standard deviation: %.3f')
		print(perr)
		# Calculate degrees of freedom of fit
		dof = len(radius_PSF) - len(p)
						   
		# Calculate best fit model
		y_fit = np.zeros(len(radius_PSF))
		for jbin in range(len(radius_PSF)):
			y_fit[jbin] = self.king_profile(radius_PSF[jbin],p[0], p[1], p[2])

		sigma_fit = round(abs(p[0]), 10)
		gamma_fit = round(abs(p[1]), 10)
		B_fit = round(abs(p[2]), 10)

		# computing the source coverage
		total_counts = np.sum(counts_PSF)
		counts_interp_norm = np.cumsum(counts_PSF/total_counts)
		interp_cum = interp1d(radius_PSF, counts_interp_norm)
		source_coverage = interp_cum(ranal)
				
		if plot_flag:
			ax.plot(radius_PSF, density_PSF_norm, '-k', linewidth=2.5, label='IRF distribution')
			ax.plot(radius_PSF, y_fit, '-r', linewidth=2.5, label='King fit')
			ax.set_xlim(0, ranal)
			ax.set_ylabel('norm. Surface Brightness [counts sr$^{-1}$]')
			ax.set_xlabel('Radial distance [deg.]')
			ax.legend(numpoints=1)
			plt.grid()
			plt.show()

		if verbose == 1:
			print('###########################################')
			print('#              SOURCE PSF2                #')
			print('###########################################')
			print('# - Energy min [MeV] = %d'% emin)
			print('# - Energy max [MeV] = %d'% emax)
			print('# - photon index [E^-gamma] = %.2f'% gindex)
			print('# - source off-axis [deg.] = %d'% theta_sel)
			print('# - source azimuthal angle [deg.] = %d'% phi_sel)
			print('# - PSF scale factor [deg.] = %.4f'% source_coverage)
			print('###########################################')
		#print('# - source coverage [%] = ')
		#print(source_coverage*100)
		#print('# - normalized PSF [deg.] = ')
		#print(radius_PSF)
		return source_coverage
