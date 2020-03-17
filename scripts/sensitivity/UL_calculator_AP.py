"""
 UL_calculator_AP.py  -  description
 ---------------------------------------------------------------------------------------------
 Computing the upper flux limit - or sensitivity - of gamma-ray missions at a given coordinate
 using aperture photometry
 ---------------------------------------------------------------------------------------------
 Usage:
 UL_calculator_AP.py <emin> <emax> <expo_on> <radius_on> <gindex> <source_theta> <source_phi> <IRF_path> <bkg_type>
 ---------------------------------------------------------------------------------------------
 Parameters:
 - emin = minimum energy in MeV
 - emax = maximum energy in MeV
 - expo_on = exposure in cm2 s
 - radius_on = maximum angle for the region in degrees
 - gindex = photon index for the power law gamma [E^-gindex]
 - source_theta = source off-axis [deg.]
 - source_phi = source azimuthal angle [deg.]
 - IRF_path = path of the IRF files 
 - bkg_type = [0 = MODEL, 1 = COUNTED]
 ---------------------------------------------------------------------------------------------
 Additional parameters:
 - if bkg_type = 0 (MODEL)
 	 - bkg_flux = flux in [cts/cm2/s/sr]
 - if bkg_type = 1 (COUNTED)
 	 - count_type = [0 = RING, 1 = EXT]
	 - expo_off = exposure in cm2 s
	 - radius_off = maximum angle for the region in degrees
	 - bkg_c = counts from the background region (in the ring or in the full external region)
 --------------------------------------------------------------------------------------------- 
 copyright            : (C) 2019 Valentina Fioretti
 email                : valentina.fioretti@inaf.it
 ----------------------------------------------
 Usage:
 UL_calculator_AP.py 
 
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.io import fits
from matplotlib import gridspec

# Import the input parameters
arg_list = sys.argv
emin = float(arg_list[1])
emax = float(arg_list[2])
expo_on = float(arg_list[3])
radius_on = float(arg_list[4])
gindex = float(arg_list[5])
source_theta = float(arg_list[6])
source_phi = float(arg_list[7])
irf_path = arg_list[8]
bkg_type = int(arg_list[9])

if bkg_type == 0:
	bkg_flux = float(arg_list[10])
	
if bkg_type == 1:
	count_type = int(arg_list[10])
	expo_off = float(arg_list[11])
	radius_off = float(arg_list[12])
	bkg_c = float(arg_list[13])	


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

# model functions

def PowerLaw(x, slope):
	return x**(-slope)

def IntPowerLaw(e1, e2, slope):
	return (((e2**(-slope + 1.))/(-slope + 1.)) - ((e1**(-slope + 1.))/(-slope + 1.)))

# select energy range
where_band = np.where((Energy_min >= emin) & (Energy_max <= emax))
energy_band = Energy[where_band]
PSF_band = PSF[where_band]
energymin_band = Energy_min[where_band]
energymax_band = Energy_max[where_band]

int_source_band = IntPowerLaw(emin, emax, gindex)
fract_source = []
for je in range(len(energy_band)):
	int_source_bin = IntPowerLaw(energymin_band[je], energymax_band[je], gindex)
	fract_source.append(int_source_bin/int_source_band)

PSF_norm = np.average(PSF_band, weights = fract_source)

print ("###########################################")
print ("#              SOURCE PSF                 #")
print ("###########################################")
print ("# - Energy min [MeV] = %d"% emin)
print ("# - Energy max [MeV] = %d"% emax)
print ("# - photon index [E^-gamma] = %.2f"% gindex)
print ("# - source off-axis [deg.] = %d"% theta_sel)
print ("# - source azimuthal angle [deg.] = %d"% phi_sel)
print ("# - normalized PSF [deg.] = %.2f"% PSF_norm)
	

#################################### 
# computing the integral sensitivity 

N_min = [1., 5., 10.]
Sign = [2., 3., 4., 5.] # Li&Ma equation 17
Omega_on = 2.*np.pi*(1. - np.cos(radius_on*(np.pi/180.)))

if bkg_type == 0:
	# parameters
	alpha = 1

	B_on = bkg_flux*Omega_on*expo_on # in counts
	N_off = B_on/alpha
	

	print ("###########################################")
	print ("#           ALGORYTHM = MODEL             #")
	print ("###########################################")
	print ("# - Exposure [cm2 s] = %.2f"% expo_on)
	print ("# - Source radius [deg] = %.2f"% radius_on)
	print ("# - Background [cts/cm2/s/sr] = %.2f"% bkg_flux)
	print ("# - alpha = %.2f"% alpha)
	print ("###########################################")

	for jmin in range(len(N_min)):

		N_source = np.float(N_min[jmin])
		for js in range(len(Sign)):

			while(1):
				if (B_on > 0.):
					N_on = N_source + alpha*N_off
					A_part = ((1. + alpha)/alpha)*(N_on/(N_on + N_off))
					B_part = (1. + alpha)*(N_off/(N_on + N_off))
					S_lima = np.sqrt(2)*math.sqrt((N_on*np.log(A_part)) + (N_off*np.log(B_part)))

					if (S_lima >= Sign[js]):
						break
					else:
						N_source += 1
				else:
					N_source = 0.
					break
	
			F_lim = (N_source/expo_on)/PSF_norm

			print ("# Sensitivity parameters:")
			print ("# - sigma = %.2f"% Sign[js])
			print ("# - N_min = %.2f"% N_min[jmin])
			print ("# Sensitivity [phot/cm2/s] = %.2e"% F_lim)
			print ("###########################################")


if bkg_type == 1:

	if count_type == 0:
		Omega_tot = 2.*np.pi*(1. - np.cos(radius_off*(np.pi/180.)))
		Omega_off = Omega_tot - Omega_on

	if count_type == 1:
		Omega_off = 2.*np.pi*(1. - np.cos(radius_off*(np.pi/180.)))

	alpha = (expo_on*Omega_on)/(expo_off*Omega_off)
	N_off = bkg_c


	print ("###########################################")
	print ("#           ALGORYTHM = COUNTED           #")
	print ("###########################################")
	print ("# - Exposure [cm2 s] = %.2f"% expo_on)
	print ("# - Source radius [deg] = %.2f"% radius_on)
	print ("# - Background [cts/cm2/s/sr] = %.2e"% N_off/expo_off/Omega_off)
	print ("# - alpha = %.2f"% alpha)
	print ("###########################################")

	for jmin in xrange(len(N_min)):

		N_source = np.float(N_min[jmin])
		for js in range(len(Sign)):

			while(1):
				if (N_off > 0.):
					N_on = N_source + alpha*N_off
					A_part = ((1. + alpha)/alpha)*(N_on/(N_on + N_off))
					B_part = (1. + alpha)*(N_off/(N_on + N_off))
					S_lima = np.sqrt(2)*math.sqrt((N_on*np.log(A_part)) + (N_off*np.log(B_part)))

					if (S_lima >= Sign[js]):
						break
					else:
						N_source += 1
				else:
					N_source = 0.
					break
	
			F_lim = N_source/expo_on

			print ("# Sensitivity parameters:")
			print ("# - sigma = %.2f"% Sign[js])
			print ("# - N_min = %.2f"% N_min[jmin])
			print ("# Sensitivity [phot/cm2/s] = %"% F_lim)
			print (F_lim)
			print ("###########################################")

