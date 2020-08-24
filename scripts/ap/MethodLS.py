# DESCRIPTION
#       Agileap: AGILE Observatory Aperture Photometry Analysis
# NOTICE
#      Any information contained in this software
#      is property of the AGILE TEAM and is strictly
#      private and confidential.
#      Copyright (C) 2005-2020 AGILE Team.
#          Bulgarelli Andrea <andrea.bulgarelli@inaf.it>
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
import matplotlib.pyplot as plt
import scipy.signal as signal
from astropy.timeseries import LombScargle

#Lomb-Scargle methods
class MethodLS:

	def __init__(self, tstartA, res, apfilename=''):
		self.tstartA = tstartA
		self.res = res
		self.apfile = apfilename
		return

	def calculateLS(self, verbose=0, plot=1, rescol=0, minfreq=-1, maxfreq=-1):
		#normalization standard model log psd
		#, normalization='standard'
		ls = LombScargle(self.tstartA, self.res[:,rescol])
		if minfreq != -1:
			#10e-6, 10e-4 , samples_per_peak=10
			frequency, power = ls.autopower(minimum_frequency=minfreq, maximum_frequency=maxfreq)
		else:
			frequency, power = ls.autopower()
		pls = ls.false_alarm_probability(power.max(), method='baluev')
		ind = power.argmax()
		maxf = frequency[ind]
		daymax = 1 / maxf / 86400.

		if verbose == 1:
			print('Max value col %2d: FAR %.5e / power %.4f / frequency %.5e Hz (%.5f days)' %(rescol, pls, power.max(), maxf, daymax) )

		if plot > 0:
			plt.subplot(1, 1, 1)
			plt.plot(frequency, power)
			plt.gca().set_xscale("log")
			plt.gca().set_yscale("log")
			#plt.subplot(3, 1, 2)
			#best_frequency = frequency[np.argmax(power)]
			#t_fit = np.linspace(0, 1)
			#y_fit = LombScargle(tstartA, ctsdataA).model(t_fit, best_frequency)
			#plt.plot(t_fit, y_fit)
		if plot == 1:
			plt.show()
		if plot == 2:
			plt.savefig(self.apfile + ".apres.ls."+str(rescol)+".png", format='png', dpi=1200)

		return pls, power.max(), maxf

	def scanLS(self, minfreq=-1, maxfreq=-1, rangemin=0, rangemax=23):
		print("freqmin %.5e"%minfreq)
		print("freqmax %.5e"%maxfreq)
		fileclean = open(self.apfile + ".apres","w")
		for i in range(rangemin,rangemax):
			pls, pmax, maxf = self.calculateLS(0 ,0, i, minfreq, maxfreq)
			daymax = 1 / maxf / 86400.
			#print("LS " + str(i) + ' ' + str(pls) + ' ' + str(pmax) + ' ' + str(maxf) + ' ' + str(daymax))
			print('Max value %2d: FAR %.5e / power %.4f / frequency %.5e Hz (%.5f days)' %(i, pls, pmax, maxf, daymax) )
			fileclean.write("LS " + str(i) + ' ' + str(pls) + ' ' + str(pmax) + ' ' + str(maxf) + ' ' + str(daymax) + "\n")
		fileclean.close()
		return

	def calculateLSexp(self, verbose=0, plot=1, rescol=3):
		ls = LombScargle(self.tstartA, self.res[:,rescol])
		frequency, power = ls.autopower()
		pls = ls.false_alarm_probability(power.max(), method='baluev')
		ind = power.argmax()
		maxf = frequency[ind]
		daymax = 1 / maxf / 86400.
		if verbose == 1:
			print('Max value '+ str(rescol) +': pls ' + str(pls) + ' with power of ' + str(power.max()) + ' at frequency ' + str(maxf) + ' Hz (' + str(daymax) + ') days' )

		if plot == 1:
			plt.subplot(1, 1, 1)
			plt.plot(frequency, power)
			plt.gca().set_xscale("log")
			plt.gca().set_yscale("log")
			#plt.subplot(3, 1, 2)
			#best_frequency = frequency[np.argmax(power)]
			#t_fit = np.linspace(0, 1)
			#y_fit = LombScargle(tstartA, ctsdataA).model(t_fit, best_frequency)
			#plt.plot(t_fit, y_fit)
			plt.show()

		return

