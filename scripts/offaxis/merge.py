import os
import pylab
import matplotlib
from matplotlib import *
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
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

    def Plotmerge(self, show=False, showcnt=True, mode="all"):

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
            ax.plot(agl_filt - self.t0, agl_sep_filt, color='blue', label='AGILE')
        if(mode=="fermi" or mode=="all"):
            ax.plot(lat_filt - self.t0, lat_sep_filt, color='red', label='Fermi-LAT')

        agilecnt_mjd = self.PlotAgileCounts()
        if showcnt==True:
            [pylab.axvline(_x, linewidth=1, color='k') for _x in agilecnt_mjd]
            for x in agilecnt_mjd:
                ax.axvline(x, linestyle='-', color='k', linewidth=1)

        ax.set_ylim(0., self.zmax+5.0)
        ax.set_xlim((self.timelimiti - self.t0)-0.2, (self.timelimitf-self.t0)+0.2)
        ax.set_xlabel('MJD')
        ax.ticklabel_format(axis="x", useOffset=False)
#        ax.set_xlabel('T - T0 [days]')
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
            f.savefig('merged_plot_'+str(self.zmax)+'_'+str(self.timelimiti)+'_'+str(self.timelimitf)+'.'+str('pdf'), format="pdf")

    def histogram_merge(self, show=False, mode="all"):

        if(mode=="agile" or mode=="all"):
            agile_center, agile_hist, agile_width = np.loadtxt('agile_histogram_visibility.txt', unpack=True)
        if(mode=="fermi" or mode=="all"):
            fermi_center, fermi_hist, fermi_width = np.loadtxt('fermi_histogram_visibility.txt', unpack=True)
        print 'Plotting histogram...'
        f = plt.figure()
        ax = f.add_subplot(111)

        if(mode=="agile" or mode=="all"):
            ax.bar(agile_center, agile_hist, align='center', color='w', edgecolor='g', width=agile_width)
        if(mode=="fermi" or mode=="all"):
            ax.bar(fermi_center, fermi_hist, align='center', color='w', edgecolor='k',width=fermi_width)

        ax.set_xlim(0., 100.)
        ax.set_ylim(0., 100.)
        ax.set_ylabel('\\% of time spent')
        ax.set_xlabel('off-axis angle $[^\\circ]$')
        labels  = [0, 10, 20, 30, 40, 50, 60, 180]
        xlabels = [0, 10, 20, 30, 40, 50, 60, 100]
        plt.xticks(xlabels, labels)
        #legend = plt.legend(loc='lower right', shadow=True, fontsize='large')

        if show==True:
            f.show()
        else:
            f.savefig('histogram_plot.png')
            f.savefig('histogram_plot.pdf', format="pdf")


    def PlotAgileCounts(self):

        agile_mjd = np.array([54732.16,54732.33]) - self.t0
        print agile_mjd
        return agile_mjd
