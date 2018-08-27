"""
    apsimlc.py
    ---------------------------------------------------------------------------------
    python script to generate a simulated light curve
    from a given .ap file with additional input from gas map
    ---------------------------------------------------------------------------------
    Dependencies:
    - python 2.7
    - numpy
    - astropy
    ---------------------------------------------------------------------------------
    Example:
    import apsimlc as aps
    calculateAPLC()
    ---------------------------------------------------------------------------------
    - background calculation based on xphs_vf.py - 2018/08/20

"""

import string, os, sys
import numpy as np
import sys
import math
from scipy.stats import norm
import matplotlib.pyplot as plt

def add_sinusoidal_signal(apfile,plot,period,phi,peak_size, delta):
    tstart = 0
    tstop = 0
    nlines = 0
    with open(apfile, "r") as ins:
       for line in ins:
          if(line != ""):
             val = line.split(" ")
             if(tstart == 0):
                tstart  = float(val[0])
             tstop   = float(val[1])
             nlines += 1
						
    tzero = tstart
    bin_number = tstop - tstart
    print('number of bins: ' + str(bin_number))
    nintervals = int(nlines)
    print('number of intervals: ' + str(nlines))
    light_curve = np.full(bin_number,0)
    light_curve_mean = np.full(bin_number,0)
    light_curve_mean_period = np.full(nintervals,0)
    omega = 2*np.pi/period
    for i in range(int(tstart),int(tstop)):
		
		lam = math.fabs(peak_size * np.sin(omega*(i) + phi) + delta)
		light_curve[i-tzero] += lam

    nline = 0
    with open(apfile, "r") as ins:
       for line in ins:
          if(line != ""):
            val = line.split(" ")
            tstart  = float(val[0])
            tstop   = float(val[1])
            meanflux = 0
            for i in range(int(tstart),int(tstop)):
               meanflux += light_curve[i-tzero]
            meanflux /= int(tstop-tstart)
            print(str(tstart) + '-' + str(tstop) + ' TT ' + str(tstop-tstart) + ' s ' + str(meanflux) + ' cts / cm2 s')
            for i in range(int(tstart),int(tstop)):
               light_curve_mean[i-tzero] = meanflux
            light_curve_mean_period[nline] = meanflux
            nline += 1

    if(plot == 1):
		#plt.plot(light_curve,'-o' ,markersize=1)
        plt.plot(light_curve_mean,'-o' ,markersize=1)
			#plt.plot(light_curve_mean_period,'-o')
        plt.savefig(apfile+".lc.png")
		#plt.show()

    return light_curve_mean_period


def calculateAPLC(help=1,  apfile='', gasvalue=-1, ranal= -1, gal = 0.7, iso = 10., emin = 100., emax = 50000.):
    if (help == "1"):
        print('--------------------------------------------------------------')
        print('                      apsimlc parameters                      ')
        print('--------------------------------------------------------------')
        print('- help: set to 1 to get usage instructions')
        print('1 apfile = file name for the AP file')
        print('2 gasvalue = the value of the gas map')
        print('3 ranal = radius of the aperture photometry')
        print('4 gal = background gal value')
        print('5 iso = background iso value (the value is implicitly multiplied by e-05)')
        print '6 emin = minimum energy [MeV]'
        print '7 emax = maximum energy [MeV]'
		#print '- gaspixelsize = the pixel size'
    else:
		
		gasvalue = float(gasvalue) #andare direttamente nei file .disp.conv.sky.gz del modello e prendere il valore da li'
        ranal = float(ranal)
        gal = float(gal)
        iso = float(iso)
        emin = float(emin)
        emax = float(emax)

		#gasdata,gashdr= fits.getdata(gasmap,header=True)
		#dx=(gashdr['CDELT1'])*(np.pi/180.) # pixel side in rad.
		#dy=(gashdr['CDELT2'])*(np.pi/180.) # pixel side in rad.

		#pix_size_sr = abs(gaspixelsize*gaspixelsize)


		#data_shape = gasdata.shape
		#ydim = data_shape[0]
		#xdim = data_shape[1]
        if emin < 1700.: psf = 0.3 # deg.
        if emin < 1000.: psf = 1. # deg.
        if emin < 400.: psf = 3.5 # deg.
        if emin < 100.: psf = 5.0 # deg. da provare
        print('selected psf: ' + str(psf))
        print('selected ranal: ' + str(ranal))
		
		#to take into account the extension of the PSF
        fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranal))
        print('fluxscalefactor: ' + str(fluxscalefactor))

		# ON - PSF region
        omega_ranal = 2.*np.pi*(1. - np.cos(ranal*(np.pi/180.)))
        e_mean = 10.**(np.log10(emin) + (np.log10(emax) - np.log10(emin))/2.)

        # Open ap file
        file = open(apfile + ".sim","w")
        file2 = open(apfile + ".sim.ap","w")
        fileclean = open(apfile + ".clean.ap","w")
        meanratetrue = 0
        meanratetrueerr = 0
        meanratesim = 0
        meanratesimerr = 0
        nlines = 0
        expsum = 0

        #clean ap file
        with open(apfile, "r") as ins:
           for line in ins:
              if(line != ""):
                 val = line.split(" ")
                 tstart  = float(val[0])
                 tstop   = float(val[1])
                 expdata = float(val[2])
                 ctsdata = int(val[3])
                 if(expdata != 0):
                    fileclean.write(str(tstart) + " " + str(tstop) + " " + str(expdata) + " " + str(ctsdata) + "\n")
        fileclean.close()
        apfile = apfile + ".clean.ap"

        #generate simulated lc
        period = 1.88 * 86400.
        phi = 0
        peak_size = 200e-08
        deltaflux = 0e-08
        lc = add_sinusoidal_signal(apfile,1, period,phi,peak_size,deltaflux)
		
        nline = 0
        with open(apfile, "r") as ins:
            for line in ins:
                if(line != ""):
                    #line_clean = line.replace("\"","").replace("\n","")
                    val = line.split(" ")
                    tstart  = float(val[0])
                    tstop   = float(val[1])
                    expdata = float(val[2]) #[cm^2 s]
					
					###############
					#expdata = 15926730
					
					#expsum = expsum + expdata
						#expdata = 5.31681e+07
                    ctsdata = float(val[3])
                    meanratetrue = meanratetrue + ctsdata / expdata
                    meanratetrueerr = meanratetrueerr + ctsdata
                    nlines = nlines + 1
                    # total background map
				    # gasvalue -> [cts] / [cm^2 s sr] #0.00061524 -> si può prendere direttamente dalld disp.conv.sky.gz 
				    # isovalue -> 10^{-5} [cts] / [cm^2 s sr]
					# gal is adimensional	
				    #bkgdata = gasdata[jy, jx]*gal + iso*(10.**(-5)) # counts/cm2/s/sr
                    bkgdata = gasvalue*gal + iso*(10.**(-5)) # cts / [cm2 s sr]
                    ctsgal = gasvalue*gal * omega_ranal*expdata #[cts]
                    ctsiso = iso*(10.**(-5)) * omega_ranal*expdata #[cts]

    				# ON background counts
    				#bkg_ON = bkgdata*omega_psf*(expdata[jy, jx]/pix_size_sr) # counts
                    bkg_ON = bkgdata*omega_ranal*expdata # cts / [cm^2 s sr] * [sr] * [cm^2 s] = cts

					#int_source_count = ((emax**(1. - index))/(1. - index) - (emin**(1. - index))/(1. - index))*(expdata[jy, jx]/pix_size_sr)*omega_psf # photons
                    index = 1.7
                    #int_source_flux = ((emax**(1. - index))/(1. - index) - (emin**(1. - index))/(1. - index))*omega_ranal # [cts / cm2 s]
					
                    fluxsource = lc[nline] * fluxscalefactor #[cts / cm2 s]
                    print(str(tstart) + '-' + str(tstop) + ' ' +  str(lc[nline]) )
                    fluxsource = 100e-08 * fluxscalefactor #[cts / cm2 s]
                    src_ON = expdata * fluxsource # [cm^2 s] * [cts] / [cm^2 s] = [cts]
                    #src_ON2 = expdata * int_source_flux # [cm^2 s] * [cts] / [cm^2 s] = [cts]

                    ctstot = bkg_ON + src_ON # [cts]
                    snr = src_ON / math.sqrt(ctstot)
                    ctsdata = np.random.poisson(ctstot)
#ctsdata = rate
                    meanratesim = meanratesim + ctsdata / expdata
                    meanratesimerr = meanratesimerr + ctsdata
                    row = str(tstart) + " " + str(tstop) + " " + str(expdata) + " " + str(ctsdata) + " " + str(bkg_ON) + " " + str(ctsgal) + " " + str(ctsiso) + " " + " " + str(src_ON) + ' ' + str(snr) + "\n"
                    file.write(row)
                    row2 = str(tstart) + " " + str(tstop) + " " + str(expdata) + " " + str(ctsdata) + "\n"
                    file2.write(row2)
                    nline += 1
	
        file.close()
        file2.close()
        print('mean rate true: ' + str(meanratetrue/nlines) + ' +/- ' + str(math.sqrt(meanratetrueerr) / expsum))
        print('mean rate sim : ' + str(meanratesim/nlines ) + ' +/- ' + str(math.sqrt(meanratesimerr) / expsum))

if __name__ == '__main__':

	# Run binned in-memory pipeline
	calculateAPLC(0, sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],100,50000)
