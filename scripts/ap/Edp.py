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
import numpy as np
import os
from astropy.io import fits

class Edp:

	def __init__(self):
		self.edpread = False
		return

	def getVal(self, trueE,  obsE,  theta,  phi):

		l = np.where(self.m_edptrueenergy==trueE)[0][0]
		m = np.where(self.m_edpobsenergy==obsE)[0][0]
		n = np.where(self.m_edptheta==theta)[0][0]
		p = np.where(self.m_edpphi==phi)[0][0]
		#8 19 16 16
		#phi, theta, obs en, true en
		return self.m_edpgrid[p,n,m,l]


	def readData(self, verbose=0):
		fits_name = os.environ['AGILE']+"/model/scientific_analysis/data/AG_GRID_G0017_SFMG_H0025.edp.gz"
		
		hdulist_edp = fits.open(fits_name)
		wcs_tab_edp = hdulist_edp[2].data
		m_edptrueenergy = wcs_tab_edp.field('TRUE_ENERGY')[0]
		m_edpobsenergy = wcs_tab_edp.field('OBS_ENERGY_CHANNEL')[0]
		m_edptheta = wcs_tab_edp.field('POLAR_ANGLE')[0]
		m_edpphi = wcs_tab_edp.field('AZIMUTH_ANGLE')[0]

		self.m_edptrueenergy = m_edptrueenergy
		self.m_edpobsenergy = m_edpobsenergy
		self.m_edptheta = m_edptheta
		self.m_edpphi = m_edpphi

		self.m_edpgrid = hdulist_edp[0].data
		
		self.edpread = True

		if verbose == 1:
			print("true energy " + str(len(m_edptrueenergy)))
			for i in range(0,len(m_edptrueenergy)):
				print(str(i+1) + " " + str(m_edptrueenergy[i]) + " check index "+ str(np.where(m_edptrueenergy==m_edptrueenergy[i])[0][0]))

			print("true m_edpobsenergy " + str(len(m_edpobsenergy)))
			for i in range(0,len(m_edpobsenergy)):
				print(str(i+1) + " " + str(m_edpobsenergy[i]) + " check index "+ str(np.where(m_edpobsenergy==m_edpobsenergy[i])[0][0]))

			print("true m_edpphi " + str(len(m_edpphi)))
			for i in range(0,len(m_edpphi)):
				print(str(i+1) + " " + str(m_edpphi[i]) + " check index "+ str(np.where(m_edpphi==m_edpphi[i])[0][0]))

			print("true m_edpphi " + str(len(m_edpphi)))
			for i in range(0,len(m_edpphi)):
				print(str(i+1) + " " + str(m_edpphi[i]) + " check index "+ str(np.where(m_edpphi==m_edpphi[i])[0][0]))




	def UpdateNormPL(self, eMin, eMax, index):

		index = 1.0-index
		m_normFactor = (pow(eMin, index)-pow(eMax, index))

		return m_normFactor

	def detCorrectionSpectraFactorSimple(self, eMin, eMax, par1, verbose=0):
	
		if(self.edpread == False):
			self.readData(verbose)
			
		expindex = 2.1

		m_edptrueenergy = self.m_edptrueenergy
		m_edpobsenergy = self.m_edpobsenergy
		m_edptheta = self.m_edptheta
		m_edpphi = self.m_edpphi

		eneChanCount = len(m_edptrueenergy)

		#print(str(m_edptrueenergy[iMin])+" "+str(m_edptrueenergy[iMax]))
		normsumpl = 0.0
		normsumple = 0.0
		iMin = 0
		iMax = 0
		for etrue in range(0,eneChanCount-1):
			if eMin == m_edptrueenergy[etrue]:
				iMin = etrue
			if eMax == m_edptrueenergy[etrue]:
				iMax = etrue
			
		if verbose == 1:
			print('EDP Energy min = %.2f'% m_edptrueenergy[iMin])
			print('EDP Energy max = %.2f'% m_edptrueenergy[iMax])

		for i in range (iMin,iMax):
			lastenergy = m_edptrueenergy[i+1]
			#print(i)
			udp1 = 0

			udp1 = self.UpdateNormPL(m_edptrueenergy[i], lastenergy, par1)
			normsumple += udp1

			udp1 = self.UpdateNormPL(m_edptrueenergy[i], lastenergy, expindex)
			normsumpl += udp1

		#print("A "+str(normsumpl)+" "+str(normsumple));
		edpArr =np.zeros(eneChanCount)

		thetaind=0
		phiind=0

		avgValuePL = 0.0
		avgValuePLE = 0.0
		for etrue in range(0,eneChanCount-1):

			lastenergy = m_edptrueenergy[etrue+1]


			for eobs in range(iMin,iMax+1):

				#print("edp: " +str(thetaind)+ " "+str(phiind)+ " "+str(etrue)+ " "+str(eobs)+ " " + str(m_edptrueenergy[etrue]) + " " + str(m_edpobsenergy[eobs]) + " " + str(m_edptheta[thetaind]) + " " + str(m_edpphi[phiind]) + " " + str(edpGrid.getVal(m_edptrueenergy[etrue], m_edpobsenergy[eobs], m_edptheta[thetaind], m_edpphi[phiind])))
				edpArr[etrue] += self.getVal(m_edptrueenergy[etrue], m_edpobsenergy[eobs], m_edptheta[thetaind], m_edpphi[phiind]) #CORRETTO


			avgValuePL  += edpArr[etrue] * self.UpdateNormPL(m_edptrueenergy[etrue], lastenergy, expindex) * 1
			avgValuePLE += edpArr[etrue] * self.UpdateNormPL(m_edptrueenergy[etrue], lastenergy, par1) * 1

		avgpl = 0.0
		avgple = 0.0


		avgpl = avgValuePL/normsumpl
		avgple = avgValuePLE/normsumple

		corr = avgpl/avgple
		#print("A " + str(m_edptheta[thetaind]) + " " + str(m_edpphi[phiind]) + " " + str(avgpl) + " " + str(avgple) + " " + str(avgpl - avgple) + " PL/" + str(avgpl / avgple))
		return corr



if __name__ == "__main__":

	par1 = 1.7
	par2 = 3403 #ec
	par3 = 0
	index = 2.1 #index, gamma1
	typefun = 0

	# reading the energy dispersion file

	#edpGrid = EdpGrid()
	

	#detCorrectionSpectraFactorSimple(4, 12, par1)
