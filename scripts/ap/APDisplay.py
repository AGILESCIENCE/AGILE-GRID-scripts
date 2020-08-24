# DESCRIPTION
#       Agileap: AGILE Observatory Aperture Photometry Analysis
# NOTICE
#      Any information contained in this software
#      is property of the AGILE TEAM and is strictly
#      private and confidential.
#      Copyright (C) 2005-2020 AGILE Team.
#          Bulgarelli Andrea <andrea.bulgarelli@inaf.it>
#          Antonio Addis <antonio.addis@inaf.it>
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


import numpy as np
from astropy.stats import bayesian_blocks
import matplotlib.pyplot as plt
import pandas as pd
#import seaborn as sb
from scipy.stats import norm


class APDisplay:

	def __init__(self, fig, axes):
		self.fig = fig
		self.axes = axes
		return
	
	#col1: histo on cts
	def plotrow(self, rown, filename, colorhist, data, xmaxcts, xmaxrate, nbinsrate):
		nbinscts = xmaxcts
		nbins = nbinsrate+1
		
		#col1: cts
		#sb.distplot(binned_2018dq1_43200["cts"], bins=nbins, kde=True, ax=self.axes[0, 0])
		n, bins, patches = self.axes[rown,0].hist(data["cts"], nbinscts, density=False, color=colorhist[0], range=[0, xmaxcts], alpha=0.5)
		binspdf = np.linspace(0, xmaxcts, xmaxcts*10.0)
		(meancts, stdcts) = norm.fit(data["cts"])
		lin = norm.pdf(binspdf, meancts, stdcts)
		self.axes[rown, 0].plot(binspdf, lin, alpha=0.9)
		
		
		n, bins, patches = self.axes[rown,0].hist(data["cts_rateWeightedMeanR4"], nbinscts, density=False, color=colorhist[1], range=[0, xmaxcts], histtype='barstacked', alpha=0.5)
		binspdf = np.linspace(0, xmaxcts, xmaxcts*10.0)
		(meanctsm, stdctsm) = norm.fit(data["cts_rateWeightedMeanR4"])
		self.axes[rown, 0].plot(binspdf, norm.pdf(binspdf, meanctsm, stdctsm), alpha=0.5)
		
		self.axes[rown, 0].set_title("cts measured (cts) ("+colorhist[0]+") / cts bkg_model (cts_rateWeightedMeanR4) ("+colorhist[1]+")")
		self.axes[rown, 0].legend(("mean: " + str(meancts.round(5)) + " std: " + str(stdcts.round(5)), "mean: " + str(meanctsm.round(5)) + " std: " + str(stdctsm.round(5))))



		#col2: rate
		n, bins, patches = self.axes[rown, 1].hist(data["rate"]*1e8, nbins, density=True, label="rate", range=[-100, xmaxrate], color=colorhist[0], alpha = 0.2, histtype='bar', linewidth=2)
		meanrate, stdrate = norm.fit(data["rate"]*1e8)
		bins2 = np.linspace(min(bins), max(bins), 100)
		self.axes[rown, 1].plot(bins2, norm.pdf(bins2, meanrate, stdrate))
		#self.axes[rown, 1].set_xlabel(r"$10^{-8} \: ph/cm^2 s$")

		
		#col2: fluxRate
		n, bins, patches = self.axes[rown, 1].hist(data["ratediffR4"]*1e8, nbins, density=True, label="rate", range=[-100, xmaxrate], color=colorhist[1], alpha = 1, histtype='step')
		meanrate2, stdrate2 = norm.fit(data["ratediffR4"]*1e8)
		bins2 = np.linspace(min(bins), max(bins), 100)
		self.axes[rown, 1].plot(bins2, norm.pdf(bins2, meanrate2, stdrate2))


	
		#col2: flux_ratediffR4
		n, bins, patches = self.axes[rown, 1].hist(data["flux_ratediffR4"]*1e8, nbins, density=True, label="flux", range=[-100, xmaxrate], color=colorhist[2], alpha=1, histtype='step')
		meanflux, stdflux = norm.fit(data["flux_ratediffR4"]*1e8)
		bins2 = np.linspace(min(bins), max(bins), 100)
		self.axes[rown, 1].plot(bins2, norm.pdf(bins2, meanflux, stdflux))
		self.axes[rown, 1].set_xlabel(r"$10^{-8} \: ph/cm^2 s$")
		
	#	self.axes[rown, 1].legend(("mean: "+'{:.2e}'.format(meanrate)+" std: "+'{:.2e}'.format(stdrate), "mean: "+'{:.2e}'.format(meanrate2)+" std: "+'{:.2e}'.format(stdrate2), "mean: "+'{:.2e}'.format(meanflux)+" std: "+'{:.2e}'.format(stdflux)))
		#self.axes[rown, 1].legend(("mean: "+'{:.2e}'.format(meanrate)+" std: "+'{:.2e}'.format(stdrate), "mean: "+'{:.2e}'.format(meanflux)+" std: "+'{:.2e}'.format(stdflux)))
		print("mean: "+'{:.2f}'.format(meanflux)+" std: "+'{:.2f}'.format(stdflux))
		#self.axes[rown, 1].legend("mean: "+'{:.2f}'.format(meanflux)+" std: "+'{:.2f}'.format(stdflux), loc="upper left")
		
		self.axes[rown, 1].set_title("rate  ("+colorhist[0]+") / rate bkg-subtracted (ratediffR4) ("+colorhist[1]+") / flux  ("+colorhist[2]+")")
		self.axes[rown, 1].legend(["mean: "+'{:.2f}'.format(meanrate)+" std: "+'{:.2f}'.format(stdrate), "mean: "+'{:.2f}'.format(meanrate2)+" std: "+'{:.2f}'.format(stdrate2), "mean: "+'{:.2f}'.format(meanflux)+" std: "+'{:.2f}'.format(stdflux)])

