#!/usr/bin/env python

import argparse
import numpy as np
import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord, FK5, Galactic
from pyrr import Quaternion, Vector3
from math import acos, atan2
from numpy.linalg import norm
import sys

def get_theta_phi(time_tt, ra, dec, logfile):

    def cartesian_to_spherical(v):
        r = norm(v)
        unit = v / r
        theta = acos(unit[2])
        phi = atan2(unit[1], unit[0])
        return theta /np.pi * 180.,  phi / np.pi * 180.

    hdulist = pyfits.open(logfile)
    tbdata = hdulist[1].data
    time = tbdata.field('TIME')
    q1 = tbdata.field('q1')
    q2 = tbdata.field('q2')
    q3 = tbdata.field('q3')
    q4 = tbdata.field('q4')

    idx = (np.abs(time - time_tt)).argmin()

    ecoord = SkyCoord(ra, dec, frame=FK5, unit="deg")
    v = Vector3([ecoord.cartesian.x, ecoord.cartesian.y, ecoord.cartesian.z])
    q = Quaternion([q1[idx], q2[idx], q3[idx], q4[idx]])
    outv = q * v
    theta, phi = cartesian_to_spherical([-outv.x, outv.z, outv.y])

    return theta, phi

def get_log(index, time):
    with open(index, "r") as f:
        for line in f:
            cols = line.split(' ')
            if float(time) >= float(cols[1]) and float(time) < float(cols[2]):
                 return cols[0]
    return None

def gal_to_eq(l, b):
    c = SkyCoord(l=l, b=b, frame=Galactic, unit="deg")
    return c.fk5.ra.deg, c.fk5.dec.deg

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "Wrong number of parameters"
        print "Usage: "+sys.argv[0]+" time_tt l b"
        print "Ex. "+sys.argv[0]+" 347470746.0 45.13 33.00"
        sys.exit(0)
    time = float(sys.argv[1])
    logfile = get_log("/AGILE_PROC3/DATA_ASDC2/INDEX/LOG.log.index", time)
    ra, dec = gal_to_eq(sys.argv[2], sys.argv[3])
    theta, phi = get_theta_phi(time, ra, dec, logfile)
    print theta, phi
