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

from math import pow
import math
import numpy as np
import os
from astropy.io import fits


class APSignificance:

	def __init__(self):
		#self.edpread = False
		return
		
	def SNR(self, verbose=0, ctsS=0, ctsB=0):
	
		snr = ctsS / math.sqrt(ctsS + 2*ctsB)
		
		return snr
	
	#alpha of Li&Ma. If not specified, alpha = ranalS / ranalB
	#ranalS = radius of analysis of S
	#ranalB = radius of analysis of B. 
	# If the background is evaluated with AGILE MLE, ranalB=10, that is the usual radius of analysis used for the evalutation of gal and iso coefficients
	def lima(self, verbose=0, alpha=-1, N_on = 0, N_off = 0, ranalS=-1, ranalB=10):

		if alpha == -1:
			alpha = ranalS / ranalB

		A_part = ((1. + alpha)/alpha)*(N_on/(N_on + N_off))
		B_part = (1. + alpha)*(N_off/(N_on + N_off))
		S_lima = np.sqrt(2)*math.sqrt((N_on*np.log(A_part)) + (N_off*np.log(B_part)))
		TS_lima = S_lima * S_lima
		if verbose == 1:
			print('######### lima ############')
			print('alpha:  %.3f' % alpha)
			print('N_OFF:  %.3f' % N_off)
			print('N_ON:   %.3f' % N_on)
			print('S_lima: %.3f' % S_lima)
		
		return S_lima
		
	def Sa(self, verbose=0, ctsTOT = 0, ctsB=0):
	
		ts = -2 * np.log(np.exp(ctsTOT-ctsB) * np.power(ctsB / ctsTOT, ctsTOT))
		s = np.sqrt(ts)
		return s


