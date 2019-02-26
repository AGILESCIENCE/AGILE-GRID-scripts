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

import scipy.signal as signal
from astropy.stats import LombScargle
import string, os, sys
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from EvalRates import *

class SimAP:

    def __init__(self):
        return

    def add_sinusoidal_signal(self, apfile,plot,period,phi,peak_size, deltaflux):
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
        bin_number = int(tstop - tstart)
        print('number of bins: ' + str(bin_number))
        nintervals = int(nlines)
        print('number of intervals: ' + str(nintervals))
        light_curve = np.zeros(bin_number)
        #light_curve_mean = np.full(int(bin_number),0)
        light_curve_mean_period = np.zeros(nintervals)
        omega = 2*np.pi/period
        print("start light_curve fill")
        for i in range(int(tstart),int(tstop)):

            lam = math.fabs(peak_size * np.sin(omega*(i) + phi) + deltaflux)

        #    lam = peak_size * np.sin(omega*(i) + phi) + deltaflux
            light_curve[int(i-tzero)] = lam
            #print(str(i) + " " + str(light_curve[int(i-tzero)]))

        print("end light_curve fill")

        nline = 0
        with open(apfile, "r") as ins:
           for line in ins:
              if(line != ""):
                val = line.split(" ")
                tstart  = float(val[0])
                tstop   = float(val[1])
                meanflux = 0
                for i in range(int(tstart),int(tstop)):
                   meanflux += light_curve[int(i-tzero)]
                meanflux /= int(tstop-tstart)
                #print(str(tstart) + '-' + str(tstop) + ' TT ' + str(tstop-tstart) + ' s ' + str(meanflux) + ' cts / cm2 s')
                #for i in range(int(tstart),int(tstop)):
                #   light_curve_mean[int(i-tzero)] = meanflux
                #print(meanflux)
                light_curve_mean_period[nline] = meanflux
                nline += 1
                #print(nline)

        if(plot == 1):
            plt.plot(light_curve_mean_period,'-o' ,markersize=1)
            #plt.plot(light_curve_mean,'-o' ,markersize=1)
            #plt.plot(light_curve_mean_period,'-o')
            plt.savefig(apfile+".lc.png")
            #plt.show()

        return light_curve_mean_period


    # RES_600_R2.ap 2 0.0006 0.7 10 100 10000
    def calculateAPLC(self, help=1, apfile='', ranal= -1, gasvalue=-1, gal = 0.7, iso = 10., emin = 100., emax = 50000., instrumentID=0):
        if (help == "1"):
            print('--------------------------------------------------------------')
            print('                      apsimlc parameters                      ')
            print('--------------------------------------------------------------')
            print('- help: set to 1 to get usage instructions')
            print('1 apfile = file name for the AP file')
            print('2 ranal = radius of the aperture photometry')
            print('3 gasvalue = the value of the gas map')
            print('4 gal = background gal value')
            print('5 iso = background iso value (the value is implicitly multiplied by e-05)')
            print('6 emin = minimum energy [MeV]')
            print('7 emax = maximum energy [MeV]')
            #print '- gaspixelsize = the pixel size'
        else:

            gasvalue = float(gasvalue) #andare direttamente nei file .disp.conv.sky.gz del modello e prendere il valore da li'
            ranal = float(ranal)
            gal = float(gal)
            iso = float(iso)
            emin = float(emin)
            emax = float(emax)

            erate = EvalRates()

            print('gasvalue %d', gasvalue)
            print('gascoeff %d', gal)
            print('isocoeff %d', iso)
            print('emin %d', emin)
            print('emax %d', emax)

            #gasdata,gashdr= fits.getdata(gasmap,header=True)
            #dx=(gashdr['CDELT1'])*(np.pi/180.) # pixel side in rad.
            #dy=(gashdr['CDELT2'])*(np.pi/180.) # pixel side in rad.

            #pix_size_sr = abs(gaspixelsize*gaspixelsize)


            #data_shape = gasdata.shape
            #ydim = data_shape[0]
            #xdim = data_shape[1]
            psf = erate.getInstrumentPSF(instrumentID, emin)

            print('selected psf: ' + str(psf))
            print('selected ranal: ' + str(ranal))

            #to take into account the extension of the PSF
            fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranal))
            print('fluxscalefactor based on PSF: ' + str(fluxscalefactor))

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
            expsum = 0

            tstartA = []
            tstopA = []
            expdataA = []
            ctsdataA = []

            #clean ap file
            with open(apfile, "r") as ins:
               for line in ins:
                  if(line != ""):
                     val = line.split(" ")

                     tstart  = float(val[0])
                     tstop   = float(val[1])
                     expdata = float(val[2])
                     ###############
                     #if(expdata != 0):
                    #     expdata = 300000

                     ctsdata = int(val[3])

                     if(expdata != 0):
                        tstartA.append(tstart)
                        tstopA.append(tstop)
                        expdataA.append(expdata)
                        ctsdataA.append(ctsdata)

                        fileclean.write(str(tstart) + " " + str(tstop) + " " + str(expdata) + " " + str(ctsdata) + "\n")
                        expsum = expsum + float(expdata)

            fileclean.close()

            apfile = apfile + ".clean.ap"

            #generate simulated lc
            period = 8.4474 * 86400.
            print('Period ',  period)
            print('Frequency ', 1 / period)
            print('Frequency2 ', 1 / (period/2))
            phi = 0
            print('Phi ', phi)
            peak_size = 50e-08
            print('peak_size ', peak_size)
            deltaflux = peak_size
            print('deltaflux ', deltaflux)


            erate.calculateRateAndSNR(0, ranal, expsum / len(expdataA), peak_size,  gasvalue, gal, iso, emin, emax, 0)

            print("Start add_sinusoidal_signal")
            lc = self.add_sinusoidal_signal(apfile, 1, period,phi,peak_size,deltaflux)
            print("End add_sinusoidal_signal")

            indexA = 0
            for x in expdataA:
                tstart  = tstartA[indexA]
                tstop   = tstopA[indexA]

                expdata = expdataA[indexA] #[cm^2 s]
                expsum = expsum + expdata

                ctsdata = ctsdataA[indexA]

                meanratetrue = meanratetrue + ctsdata / expdata
                meanratetrueerr = meanratetrueerr + ctsdata

                # total background map
                # gasvalue -> [cts] / [cm^2 s sr] #0.00061524 -> si puo' prendere direttamente dalld disp.conv.sky.gz
    		    # isovalue -> 10^{-5} [cts] / [cm^2 s sr]
    			# gal is adimensional
    		    #bkgdata = gasdata[jy, jx]*gal + iso*(10.**(-5)) # counts/cm2/s/sr
                bkgdata = gasvalue*gal + iso*(10.**(-5)) # cts / [cm2 s sr]
                ctsgal = gasvalue * gal * omega_ranal *expdata #[cts]
                ctsiso = iso*(10.**(-5)) * omega_ranal*expdata #[cts]

    			# ON background counts
    			#bkg_ON = bkgdata*omega_psf*(expdata[jy, jx]/pix_size_sr) # counts
                bkg_ON = bkgdata*omega_ranal*expdata # cts / [cm^2 s sr] * [sr] * [cm^2 s] = cts

    			#int_source_count = ((emax**(1. - index))/(1. - index) - (emin**(1. - index))/(1. - index))*(expdata[jy, jx]/pix_size_sr)*omega_psf # photons
                #index = 1.7
                #int_source_flux = ((emax**(1. - index))/(1. - index) - (emin**(1. - index))/(1. - index))*omega_ranal # [cts / cm2 s]

                #A) flux source from variabile LC
                #print(lc[indexA])
                fluxsource = lc[indexA] * fluxscalefactor #[cts / cm2 s]
                #print(str(tstart) + '-' + str(tstop) + ' ' +  str(lc[indexA]) )

                #B) flux source constant
                #fluxsource = 100e-08 * fluxscalefactor #[cts / cm2 s]

                #Calculation of counts
                src_ON = expdata * fluxsource # [cm^2 s] * [cts] / [cm^2 s] = [cts]
                #src_ON2 = expdata * int_source_flux # [cm^2 s] * [cts] / [cm^2 s] = [cts]

                ctstot = bkg_ON + src_ON # [cts]
                snr = src_ON / math.sqrt(float(ctstot))
                ctsdata = np.random.poisson(ctstot)
                #ctsdata = ctstot
                meanratesim = meanratesim + ctsdata / expdata
                meanratesimerr = meanratesimerr + ctsdata

                #decide if write real or simulated data
                writerealdata = 0
                if writerealdata == 1:
                    ctsdata = ctsdataA[indexA]

                row = str(tstart) + " " + str(tstop) + " " + str(expdata) + " " + str(ctsdata) + " " + str(fluxsource) + " " + str(bkg_ON) + " " + str(ctsgal) + " " + str(ctsiso) + " " + " " + str(src_ON) + ' ' + str(snr) + "\n"
                file.write(row)
                row2 = str(tstart) + " " + str(tstop) + " " + str(expdata) + " " + str(ctsdata) + "\n"
                file2.write(row2)

                indexA = indexA + 1

            file.close()
            file2.close()
            print('nlines ' + str(len(expdataA)))
            print('expsum ' + str(expsum))
            print('mean rate true: ' + str(meanratetrue/len(tstartA)) + ' +/- ' + str(math.sqrt(meanratetrueerr) / expsum))

            print('mean rate sim : ' + str(meanratesim/len(tstartA) ) + ' +/- ' + str(math.sqrt(meanratesimerr) / expsum))
