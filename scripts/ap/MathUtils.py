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

import numpy as np

class MathUtils:
	def steradiansCone(self, ranal):
		return 2.*np.pi*(1. - np.cos(ranal*(np.pi/180.))) #[sr]
		#VF
		#Omega_on = 2.*np.pi*(1. - np.cos(radius_on*(np.pi/180.)))
		
	def steradiansSquarePixel(self, dimpixel, theta=30):
		pixel = (np.pi/180. * dimpixel)**2
		return pixel * np.sin(np.pi/180. * theta) / (np.pi/180. * theta) #[sr]
		#omega_ranalAC = (np.pi/180. * ranal)**2 * np.sin(np.pi/180. * 30) / (np.pi/180. * 30) #[sr]
		
	def lowCountsError(self, n_i):
			if n_i + 0.25 > 0:
				s_ip = 0.5 + np.sqrt(n_i + 0.25)
				s_im = -0.5 + np.sqrt(n_i + 0.25)
			else:
				s_ip = 0
				s_im = 0
			s_irms = np.sqrt((s_ip*s_ip + s_im*s_im) / 2.0)
			return s_irms


