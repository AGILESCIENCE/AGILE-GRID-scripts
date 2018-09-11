import scipy.signal as signal
from astropy.stats import LombScargle
import string, os, sys
import numpy as np
import sys
import math
import matplotlib.pyplot as plt
from EvalRates import *

class GammaAP:

    def __init__(self):
        return

    def loadDataAPAGILE(self, apfile):
        self.diml = 0
        with open(apfile, "r") as ins:
            for line in ins:
                if(line != ""):
                    val = line.split(" ")
                    expdata = float(val[2])
                    if(expdata != 0):
                        self.diml += 1

        self.tstartA = np.zeros(self.diml)
        self.tstopA = np.zeros(self.diml)
        self.expdataA = np.zeros(self.diml)
        self.ctsdataA = np.zeros(self.diml)

        nline = 0
        with open(apfile, "r") as ins:
            for line in ins:
                if(line != ""):
                    val = line.split(" ")
                    tstart  = float(val[0])
                    tstop   = float(val[1])
                    expdata = float(val[2])
                    ctsdata = float(val[3]) #select the right column
                    if(expdata != 0):
                        self.tstartA[nline] = tstart
                        self.tstopA[nline] = tstop
                        self.expdataA[nline] = expdata
                        self.ctsdataA[nline] = ctsdata
                        nline += 1

    def normalizeAP(self, apfile):
        return;

    def calculateLS(self, verbose=0, plot=1):
        #normalization standard model log psd
        ls = LombScargle(self.tstartA, self.ctsdataA, normalization='standard')
        frequency, power = ls.autopower(minimum_frequency=1e-7, maximum_frequency=1e-4, samples_per_peak=100)
        pls = ls.false_alarm_probability(power.max(), method='baluev')
        ind = power.argmax()
        maxf = frequency[ind]
        daymax = 1 / maxf / 86400.

        if verbose == 1:
            print('pls ' + str(pls) + ' with power of ' + str(power.max()) + ' at frequency ' + str(maxf) + ' Hz (' + str(daymax) + ') days' )

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

        return pls, maxf
