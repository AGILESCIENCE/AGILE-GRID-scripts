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

"""
	apsimlc.py
	---------------------------------------------------------------------------------
	python script to generate a simulated light curve
	from a given .ap file with additional input from gas map
	---------------------------------------------------------------------------------
	Dependencies:
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
plt.switch_backend('agg')
from EvalRates import *

class SimAP:

	def __init__(self):
		return

	def find_line(self,x1,y1,x2,y2):
		m = float(y2 - y1)/(x2-x1)
		q = y1 - (m*x1)
		def line(x):
			return m * x + q
		return line

	def get_flux_from_window_file(self,time,windowfile):

		flux = -1
		with open(windowfile, "r") as file:
			for cnt, line in enumerate(file):

				if(line != ""):
					val = line.split(" ")
					tstart = float(val[0])
					tstop = float(val[1])
					minflux = float(val[2])
					maxflux = float(val[3])


					# if bin is inside the window calculate flux
					if(time>=tstart and time<=tstop):
						line = self.find_line(tstart,minflux,tstop,maxflux)
						flux = line(time)

						break

		return flux

	def add_lc_signal(self,apfile,windowfile,plot, fluxscalefactor):
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
		print("tstart="+str(tstart))
		print("tstop="+str(tstop))

		print('number of bins: ' + str(bin_number))
		nintervals = int(nlines)
		print('number of intervals: ' + str(nintervals))
		light_curve = np.zeros(bin_number)
		#light_curve_mean = np.full(int(bin_number),0)
		light_curve_mean_period = np.zeros(nintervals)


		for i in range(int(tstart),int(tstop)):

			#print("flux"+str(self.get_flux_from_window_file(i,windowfile)))

			light_curve[int(i-tzero)] = self.get_flux_from_window_file(i,windowfile) * fluxscalefactor

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
						 #print("calculate"+str(i))
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
			#plt.plot(light_curve,'-o' ,markersize=1)
			#plt.plot(light_curve_mean,'-o' ,markersize=1)
			plt.plot(light_curve_mean_period,'-o')
			plt.savefig(apfile+".lc.png")
			#plt.show()

		return light_curve_mean_period

	def add_sinusoidal_signal(self, apfile,plot,period,phi,peak_size, deltaflux, fluxscalefactor):
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

		if(period != 0):
			omega = 2*np.pi/period
		print("start light_curve fill")
		for i in range(int(tstart),int(tstop)):

			#constant model
			if(period == 0):
				lam = math.fabs(peak_size*fluxscalefactor)
			else:
				lam = math.fabs(peak_size *fluxscalefactor * np.sin(omega*(i) + phi) + deltaflux * fluxscalefactor)

			#lam = peak_size * np.sin(omega*(i) + phi) + deltaflux
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
	def calculateAPLC(self, help=1, mode='0', windowfile='', apfile='', ranal= 1, gasvalue=0.0006, gal = 0.7, iso = 10., emin = 100., emax = 10000., period = 729855.36, expdatafixed=0, peak_size = 50e-08, gindex=2.1, instrumentID=0):
		if (help == "1"):
			print('--------------------------------------------------------------')
			print('					  apsimlc parameters					  ')
			print('--------------------------------------------------------------')
			print('- help: set to 1 to get usage instructions')
			print('1 mode = 0 to simulate sinusoidal signal, 1 to simulate signal from windowfile')
			print('2 windowfile = file name for the window file')
			print('3 apfile = file name for the AP file')
			print('4 ranal = radius of the aperture photometry')
			print('5 gasvalue = the value of the gas map. NB: andare direttamente nei file .disp.conv.sky.gz del modello e prendere il valore')
			print('6 gal = background gal value')
			print('7 iso = background iso value (the value is implicitly multiplied by e-05)')
			print('8 emin = minimum energy [MeV]')
			print('9 emax = maximum energy [MeV]')
			print('10 period in seconds 8.4474 * 86400 = 729855.36. NOTE: if period=0 similation of a constant model with peak_size value')
			print('11 expdatafixed = fix exposure data to a fixed value if > 0. The value of the exposure is expdatafixed')
			print('12 peak_size = the maximum value of the sinusoidal source to be simulated. 50e-08')
			print('13 gindex = spectral index of the gamma-ray source. Default 2.1')
			#print '- gaspixelsize = the pixel size'
		else:

			gasvalue = float(gasvalue) #andare direttamente nei file .disp.conv.sky.gz del modello e prendere il valore da li'
			ranal = float(ranal)
			gal = float(gal)
			iso = float(iso)
			emin = float(emin)
			emax = float(emax)

			erate = EvalRates()

			print('gasvalue %.5f'%gasvalue)
			print('gascoeff %.2f'% gal)
			print('isocoeff %.2f'% iso)
			print('emin [MeV]: %d'% emin)
			print('emax [MeV]: %d'% emax)
			print('gindex %.2f'% gindex)
			#gasdata,gashdr= fits.getdata(gasmap,header=True)
			#dx=(gashdr['CDELT1'])*(np.pi/180.) # pixel side in rad.
			#dy=(gashdr['CDELT2'])*(np.pi/180.) # pixel side in rad.

			#pix_size_sr = abs(gaspixelsize*gaspixelsize)


			#data_shape = gasdata.shape
			#ydim = data_shape[0]
			#xdim = data_shape[1]
			#psf = erate.getInstrumentPSF(instrumentID, emin)

			#print('selected psf   [deg]: ' + str(psf))
			#print('selected ranal [deg]: ' + str(ranal))

			#to take into account the extension of the PSF
			#fluxscalefactor = math.fabs(1-2*norm(0,  psf).cdf(ranal))
			#print('fluxscalefactor based on PSF: %.5f' % fluxscalefactor)
			fluxscalefactor=erate.getFluxScaleFactor(verbose=1, gindex=gindex, ranal=ranal, emin=emin, emax=emax)
			print("fluxscalefactor= " + str(fluxscalefactor))

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
			ctsdataFinal = []

			#clean ap file
			with open(apfile, "r") as ins:
				for line in ins:
					if(line != ""):
						val = line.split(" ")

						tstart  = float(val[0])
						tstop   = float(val[1])
						expdata = float(val[2])
						###############
						if(expdatafixed > 0):
							if(expdata != 0):
								expdata = expdatafixed

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

			print("--- Signal model")
			if(period > 0):
				print('Period     [s]: %.3f '% period)
				print('Frequency [Hz]: %.3e '% (1 / period))
				#print('Frequency2 ', 1 / (period/2))
			else:
				print("Constant")
			phi = 0
			print('Phi ', phi)
			print('peak_size [ph/cm2 s]: ', peak_size)
			deltaflux = peak_size
			print('deltaflux [ph/cm2 s]: ', deltaflux)
			print("--- End signal model")


			erate.calculateRateAndSig(verbose=1, ranalS=ranal, exposure=expsum / len(expdataA), fluxsource=peak_size,  gasvalue=gasvalue, gal=gal, iso=iso, emin=emin, emax=emax, gindex=gindex, instrumentID=0)

			print("Start add signal")
			if(mode == '0'):
				#it is necessary to multiply for fluxcorrection factor
				lc = self.add_sinusoidal_signal(apfile, 1, period,phi,peak_size,deltaflux, fluxscalefactor)
			elif(mode == '1'):
				lc = self.add_lc_signal(apfile, windowfile,1, fluxscalefactor)
			print("End add signal")

			indexA = 0
			expsum = 0
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
				fluxsource = lc[indexA] * fluxscalefactor # [cts / cm2 s]
				#print(str(tstart) + '-' + str(tstop) + ' ' +  str(lc[indexA]) )

				#B) flux source constant
				#fluxsource = 100e-08 * fluxscalefactor #[cts / cm2 s]

				#Calculation of counts
				src_ON = expdata * fluxsource # [cm^2 s] * ( [cts] / [cm^2 s] ) = [cts]
				#src_ON2 = expdata * int_source_flux # [cm^2 s] * [cts] / [cm^2 s] = [cts]

				ctstot = float(bkg_ON + src_ON) # [cts]
				snr = src_ON / math.sqrt(ctstot)
				ctsdata = np.random.poisson(ctstot)

				ctsdataFinal.append(ctsdata)
				#ctsdata = ctstot
				#print(ctstot)
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

			plot = 1
			if(plot == 1):
				plt.clf()
				plt.plot(ctsdataFinal,'-o' ,markersize=1)

				plt.savefig(apfile+".lc.final.png")
				#plt.show()

			print('nlines: ' + str(len(expdataA)))
			print('expsum [cm2 s]: ' + str(expsum))
			print('mean expsum [cm2 s]: %.2f' % (expsum / indexA))
			print('mean rate true [cts / cm2 s]: %.3e +/- %.3e ' % (meanratetrue/len(tstartA), math.sqrt(meanratetrueerr) / expsum))

			print('mean rate sim  [cts / cm2 s]: %.3e +/- %.3e ' % (meanratesim/len(tstartA), math.sqrt(meanratesimerr) / expsum))


if __name__ == "__main__":

	#execute this class from command nline

	apfile = sys.argv[1]
	windowfile = sys.argv[2]
	mode = sys.argv[3]
	ranal= sys.argv[4]

	simap = SimAP()
	#simap.calculateAPLC(mode=mode, windowfile=windowfile, apfile=apfile)
	simap.calculateAPLC(mode=mode,apfile=apfile, windowfile=windowfile, ranal=ranal, peak_size=1000e-08, period=0, gindex=1.8, gasvalue=0.0006, gal = 0.7, iso = 10., emin = 100., emax = 10000.)


	#test
