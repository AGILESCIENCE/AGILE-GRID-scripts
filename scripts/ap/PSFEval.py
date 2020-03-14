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
import sys
import math
from astropy.io import fits

class PSFEval:

	def __init__(self):
		return

	# model functions
	def PowerLaw(self, x, slope):
		return x**(-slope)

	def IntPowerLaw(self, e1, e2, slope):
		return (((e2**(-slope + 1.))/(-slope + 1.)) - ((e1**(-slope + 1.))/(-slope + 1.)))


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
