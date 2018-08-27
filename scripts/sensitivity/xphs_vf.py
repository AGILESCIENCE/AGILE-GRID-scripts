"""
    xphs_vf.py
    ---------------------------------------------------------------------------------
    python script to create an upper limit flux map 
    from a given exposure and background flux
    ---------------------------------------------------------------------------------
    Dependencies:
    - python 2.7
    - numpy
    - astropy
    ---------------------------------------------------------------------------------
    Example:
    import xphs_vf as xp
    xp.computeUL()
    xp.script()
	
	- CAT2
	import xphs_vf as xs
	xs.computeUL(help = 0, sens_type = 'INT', expmap='/Users/fioretti/pro_python/Maps/POINT01C.exp', gasmap='/Users/fioretti/pro_python/Maps/POINT01C.gas', gal = 0.7, iso = 12., cont = 10., sign = 4., emin = 100., emax = 10000., model = 0, index = 2.1, par2 = 0, par3 = 0)

    ---------------------------------------------------------------------------------
    - first release by A. Giuliani (INAF IASF Milano)
    - upgrade by V. Fioretti (INAF OAS Bologna) - 2018/04/23
    
"""

import string, os, sys
from astropy.io import fits
import numpy as np
import sys

def computeUL(help=1, sens_type = 'DIFF', expmap='', gasmap='', gal = 0.7, iso = 10., cont = 10., sign = 5., emin = 100., emax = 50000., model = 0, index = 0, par2 = 0, par3 = 0):
    
    if (help == 1):
    
        print '--------------------------------------------------------------'
        print '                      xphs_vf parameters                      '
        print '--------------------------------------------------------------'
    
        print '- help: set to 1 to get usage instructions'
        print '- sens_type:'
        print '     - DIFF = differential sensitivity'
        print '     - INT = integral sensitivity'
        print '- expmap = file name for the EXP map'
        print '- gasmap = file name for the GAS map'
        print '- gal = background gal value'
        print '- iso = background iso value (the value is implicitly multiplied by e-05)'
        print '- cont = minimum number of signal photons'
        print '- sign = minimum significance'
        print '- emin = minimum energy [MeV]'
        print '- emax = maximum energy [MeV]'
        print '- model = source spectral model (0 = power-law in the form E^-index)'
        print '- index = power-law index'
        print '- par2 = model parameter 2'
        print '- par3 = model parameter 3'
    
    else:

        expdata,exphdr= fits.getdata(expmap,header=True)
        gasdata,gashdr= fits.getdata(gasmap,header=True)
        
        if emin < 400.: psf = 4. # deg.
        if emin < 1000.: psf = 2 # deg.
        if emin < 1700.: psf = 0.7 # deg.

        # ON - PSF region
        omega_psf = 2.*np.pi*(1. - np.cos(psf*(np.pi/180.)))

        # ON + OFF - ROI = 10 deg.
        roi = 10.  # deg.
        omega_roi = 2.*np.pi*(1. - np.cos(roi*(np.pi/180.)))
        
        e_mean = 10.**(np.log10(emin) + (np.log10(emax) - np.log10(emin))/2.)
        
        data_shape = expdata.shape
        ydim = data_shape[0]
        xdim = data_shape[1]

        dx=(exphdr['CDELT1'])*(np.pi/180.) # pixel side in rad.
        dy=(exphdr['CDELT2'])*(np.pi/180.) # pixel side in rad.

        pix_size_sr = abs(dx*dy)

        
        ul_data = np.zeros(shape=(ydim,xdim))
        
        for jy in xrange(ydim):
            for jx in xrange(xdim):
            	if (expdata[jy, jx] != 0):
					# source model
					if (model == 0):
						# power-law
						int_source_count = ((emax**(1. - index))/(1. - index) - (emin**(1. - index))/(1. - index))*(expdata[jy, jx]/pix_size_sr)*omega_psf # photons
						int_source_flux = ((emax**(1. - index))/(1. - index) - (emin**(1. - index))/(1. - index))*omega_psf # photons/cm2/s
						if (sens_type == 'DIFF'): source_flux = (e_mean**(-index))*omega_psf # photons/cm2/s

					# total background map
					bkgdata = gasdata[jy, jx]*gal + iso*(10.**(-5)) # counts/cm2/s/sr

					# ON background counts
					bkg_ON = bkgdata*omega_psf*(expdata[jy, jx]/pix_size_sr) # counts

					# ON background counts
					bkg_TOTAL = bkgdata*omega_roi*(expdata[jy, jx]/pix_size_sr) # counts

					# OFF background counts
					bkg_OFF = bkg_TOTAL - bkg_ON # counts

					alpha = bkg_ON/bkg_OFF
				
					#print 'alpha', alpha

					N_source = cont
					while(1):
						if (bkg_ON > 0.):
							N_on = N_source + bkg_ON
							A_part = ((1. + alpha)/alpha)*(N_on/(N_on + bkg_OFF))
							B_part = (1. + alpha)*(bkg_OFF/(N_on + bkg_OFF))
							S_lima = np.sqrt(2.)*np.sqrt((N_on*np.log(A_part)) + (bkg_OFF*np.log(B_part)))
							flux_ratio = N_source/bkg_ON
							#print 'N_source, S_lima ', N_source, S_lima
							if (S_lima >= sign):
								break
							else:
								N_source += 1
				
						else:
							N_source = 0.
							break

					if (sens_type == 'DIFF'):
						F_lim = (N_source/int_source_count)*source_flux
					if (sens_type == 'INT'):
						F_lim = (N_source/int_source_count)*int_source_flux
                else:
					F_lim = 0.
					
                ul_data[jy, jx] = F_lim # phot/cm2/s

        ulfile=expmap+'.ul'
        fits.writeto(ulfile,ul_data,exphdr, overwrite=1)




def script(string="EXP", list='list_name', sens_type = 'DIFF', expmap='', gasmap='', gal = 0.7, iso = 10., cont = 10., Sign = 5., emin = 100., emax = 50000., model = 0, index = 0, par2 = 0, par3 = 0):

    os.system("ls *"+string+" > "+list)
    reading=open(list,"r")

    for i in range(100):
        file = reading.readline()
        file = string.replace( file, '\n', '' )

    if file != "":
        computeUL(help = 0, sens_type = 'DIFF', expmap='', gasmap='', gal = 0.7, iso = 10., cont = 10., sign = 5., emin = 100., emax = 50000., model = 0, index = 0, par2 = 0, par3 = 0)
    
    os.system("rm "+list)


