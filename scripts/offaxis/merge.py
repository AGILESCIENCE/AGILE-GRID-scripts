import os
import pylab
import matplotlib
from matplotlib import *
import matplotlib.pyplot as plt
import numpy as np
import astropy
from astropy.time import Time
from astropy import coordinates as coord
from astropy import units as u
from astropy.coordinates import SkyCoord
import pyfits
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import bisect


class merge:
    """ This class provides plots to check the off-axis angle of a source wrt AGILE center FoV.    """

    def __init__(self, zmax=75., timelimiti=-1, timelimitf=-1, step=0.1, t0=0.):
        self.zmax        = zmax
        self.timelimiti  = timelimiti
        self.timelimitf  = timelimitf
        self.step        = step*10.
        self.t0          = t0

    def Plotmerge(self, show=False, showcnt=False, mode="all"):

        if(mode=="agile" or mode=="all"):
            agl_meantime, agl_separation = np.loadtxt('time_vs_separation_agile.txt', unpack=True)
        if(mode=="fermi" or mode=="all"):
            lat_meantime, lat_separation = np.loadtxt('time_vs_separation_fermi.txt', unpack=True)

        if(mode=="agile" or mode=="all"):
            agl_filt = agl_meantime[(agl_meantime > self.timelimiti) & (agl_meantime < self.timelimitf)]
            agl_sep_filt = agl_separation[(agl_meantime > self.timelimiti) & (agl_meantime < self.timelimitf)]

        if(mode=="fermi" or mode=="all"):
            lat_filt = lat_meantime[(lat_meantime > self.timelimiti) & (lat_meantime < self.timelimitf)]
            lat_sep_filt = lat_separation[(lat_meantime > self.timelimiti) & (lat_meantime < self.timelimitf)]

        print 'Plotting figure...'
        f = plt.figure()
        ax = f.add_subplot(111)
#ax.plot(agl_meantime, agl_separation, '-r', lat_meantime, lat_separation, 'gs')
#        ax.plot(agl_filt - self.t0, agl_sep_filt, '-r', label='AGILE')
#        ax.plot(lat_filt - self.t0, lat_sep_filt, 'gs', label='Fermi-LAT')
        if(mode=="agile" or mode=="all"):
            ax.plot(agl_filt - self.t0, agl_sep_filt, color='gray', label='AGILE')
        if(mode=="fermi" or mode=="all"):
            ax.plot(lat_filt - self.t0, lat_sep_filt, '+', color='red', markersize=2, label='Fermi-LAT')

##        agilecnt_mjd = self.PlotAgileCounts()
##        if showcnt==True:
##            [pylab.axvline(_x, linewidth=3, color='k') for _x in agilecnt_mjd]
#            ax.axvline(agilecnt_mjd, linestyle='-', color='k', linewidth='3')

        ax.set_ylim(0., self.zmax+5.0)
        ax.set_xlim((self.timelimiti - self.t0)-0.2, (self.timelimitf-self.t0)+0.2)
#        ax.set_xlabel('MJD')
        ax.set_xlabel('T - T0 [days]')
        ax.set_ylabel('off-axis angle [$^{\\circ}$]')

        legend = plt.legend(loc='lower right', shadow=True, fontsize='large')

        print 'Saving figure...'
#        ax.set_xlim(np.min(agl_filt)-self.t0, np.max(agl_filt)-self.t0)
        ax.set_xlim(np.min(agl_filt-self.t0), np.max(agl_filt-self.t0))
        ax.set_title(str(self.zmax)+'_'+str(self.timelimiti)+'_'+str(self.timelimitf))
#        ax.set_xlim(self.timelimiti, self.timelimitf)
#        ax.set_xlim(self.timelimiti - self.t0, self.timelimitf - self.t0)
        if show==True:
            f.show()
        else:
            f.savefig('merged_plot_'+str(self.zmax)+'_'+str(self.timelimiti)+'_'+str(self.timelimitf)+'.'+str('png'))


    def PlotAgileCounts(self):

        agile_mjd = np.array([57597.530897,57598.274395,57598.402544,57598.597396,57598.834219,57598.849626,
                              57598.924719,57599.030639,57599.108716]) - self.t0
        print agile_mjd
        return agile_mjd
