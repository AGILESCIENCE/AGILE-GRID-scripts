#!/usr/bin/env python
"""
  Copyright (c) prior to 2016 S. Cutini, A. Giuliani
  Copyright (c) 2016 S. Cutini, A. Giuliani, V. Fioretti, A. Zoli
  Copyright (c) 2017 S. Cutini, A. Giuliani, V. Fioretti, A. Zoli, A. Bulgarelli

  Dependencies:
    - python 2.7
    - numpy
    - astropy
"""

import os
import tempfile
from sys import argv, exit
from math import pow, sqrt
from numpy import zeros, cos, sin, arccos, deg2rad, rad2deg, log
import pyfits
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5

log_index = "/ASDC_PROC2/DATA_2/INDEX/LOG.log.index"
evt_index = "/ASDC_PROC2/FM3.119_2/INDEX/EVT.index"

def fromGaltoRaDec(l_coord, b_coord):
    c = SkyCoord("galactic", l=l_coord, b=b_coord, unit="deg")
    c2 = c.transform_to(FK5(equinox='J2000'))
    return c2.ra.deg, c2.dec.deg

def grb_pipe(ofd, sigmafd, evt_file='', log_file='', GRB_time = 0., GRB_l = 0., GRB_b = 0., GRB_ra=0., GRB_dec=0., t1s=0., t2s=0., t1b=0., t2b=0., raggio=0., fov=0., ea_th=0.):
    if (evt_file=='' or log_file==''):
        print 'FATAL ERROR: Wrong input parameters!!!!'
        return

    print >> ofd, ''
    print >> ofd, 'GRB T0:', GRB_time
    print >> ofd, 'GRB T0 - t1s:', GRB_time + t1s
    print >> ofd, 'GRB T1 + t2s:', GRB_time + t2s
    print >> ofd, 'GRB (L, B):', GRB_l, GRB_b
    print >> ofd, 'GRB (Ra, Dec):', GRB_ra, GRB_dec
    print >> ofd, 'Ricerca eventi (Raggio):', raggio
    print >> ofd, 'Ricerca eventi (Tmin, Tmax):', t1s, t2s
    print >> ofd, 'Background (Tmin, Tmax):', t1b, t2b
    print >> ofd, 'F.O.V.:', fov
    print >> ofd, 'Albedo cut:', ea_th
    print >> ofd, ''
    print >> sigmafd, GRB_l, GRB_b,

    # reading log file
    print "logfile: " + log_file
    hdulist_log = pyfits.open(log_file)
    tbdata_log = hdulist_log[1].data
    TIME_log = tbdata_log.field('TIME')  # sec
    RA_earth = tbdata_log.field('EARTH_RA') # deg.
    DEC_earth = tbdata_log.field('EARTH_DEC') # deg.
    livetime = tbdata_log.field('LIVETIME') # ms
    ra_punt = tbdata_log.field('ATTITUDE_RA_Y') # deg.
    dec_punt = tbdata_log.field('ATTITUDE_DEC_Y') # deg.
    phase_log = tbdata_log.field('PHASE')

    TIMEnew_log = TIME_log - GRB_time # time starting point is the source T0.

    # reading evt file
    print "evtfile: "+evt_file
    hdulist = pyfits.open(evt_file)
    tbdata = hdulist[1].data
    RAcol = tbdata.field('RA') # deg.
    DECcol = tbdata.field('DEC') # deg.
    TIMEcol = tbdata.field('TIME') # sec.
    PH_col = tbdata.field('PH_EARTH') # deg.
    EVcol = tbdata.field('EVSTATUS') # Event Classification Flag
    Enecol = tbdata.field('ENERGY') # MeV
    THETAcol = tbdata.field('THETA') # deg., coordinates in the P/L Ref Sys.
    PHASEcol = tbdata.field('PHASE') # deg.

    TIMEnew = TIMEcol - GRB_time # time starting point is the source T0.

    ra = deg2rad(RAcol)
    dec = deg2rad(DECcol)
    ra_ea = deg2rad(RA_earth)
    dec_ea = deg2rad(DEC_earth)

    Vx_ea = zeros(len(dec_ea),float)
    Vy_ea = zeros(len(dec_ea),float)
    Vz_ea = zeros(len(dec_ea),float)

    argx_ea = zeros(len(dec_ea),float)
    argy_ea = zeros(len(dec_ea),float)
    argz_ea = zeros(len(dec_ea),float)
    arg_ea = zeros(len(dec_ea),float)

    DELTA_ea = zeros(len(dec_ea),float)

    expo = zeros(len(dec_ea),float)
    expo1 = zeros(len(dec_ea),float)
    expo2 = zeros(len(dec_ea),float)
    exposure = zeros(len(dec_ea),float)

    t1b = max(min(TIMEnew),t1b)
    t2b = min(max(TIMEnew),t2b)

    t1b = max(min(TIMEnew_log),t1b)   # changing the range of the background time if outside the evt file time
    t2b = min(max(TIMEnew_log),t2b)   # changing the range of the background time if outside the evt file time

    if t2s > t2b: # adjust t2s to be <= t2b (this avoids the tback=0 and division by 0 in the following code)
        t2s = t2b
    if t1s < t1b:
        t1s = t1b

    print >> ofd, 'Nuova Ricerca eventi (Tmin, Tmax):', t1s, t2s
    print >> ofd, 'Nuovo Background (Tmin, Tmax):', t1b, t2b

    # GRB in coord. cartesiane
    Vxgrb = cos(deg2rad(GRB_dec)) * cos(deg2rad(GRB_ra))
    Vygrb = cos(deg2rad(GRB_dec)) * sin(deg2rad(GRB_ra))
    Vzgrb = sin(deg2rad(GRB_dec))

    # Distanza angolare grb-eventi
    Vx = cos(dec) * cos(ra)
    Vy = cos(dec) * sin(ra)
    Vz = sin(dec)

    argx = Vxgrb * Vx
    argy = Vygrb * Vy
    argz = Vzgrb * Vz
    arg = argx + argy + argz

    DELTA = rad2deg(arccos(arg))

    # Off-axis angle del GRB (with respect to the telescope attitude)
    Vx_punt = cos(deg2rad(dec_punt)) * cos(deg2rad(ra_punt))
    Vy_punt = cos(deg2rad(dec_punt)) * sin(deg2rad(ra_punt))
    Vz_punt = sin(deg2rad(dec_punt))

    argx_punt = Vxgrb * Vx_punt
    argy_punt = Vygrb * Vy_punt
    argz_punt = Vzgrb * Vz_punt
    arg_punt = argx_punt + argy_punt + argz_punt

    OFF = rad2deg(arccos(arg_punt))

    # Calcolo segnale e background
    source = 0
    bkg = 0

    for k in range(len(dec)):
#        if EVcol[k] != 'L':
	  if DELTA[k] <raggio:  # if the event is within the ring
            if Enecol[k] > 0.:        # if the event energy is > 0
             if PH_col[k] >ea_th:        # removing the events from albedo
                if PHASEcol[k] != 1:      # requiring PHASE not equal 1     ??????????
                  if THETAcol[k]<fov:       # theta of the event in the P/L Sys. Ref. within FOVRADMAX
                    if TIMEnew[k] > t1s:
                        if TIMEnew[k] < t2s: # if the event time is within the T1-T2 of the source
                            source = source + 1   # counting the events that pass the selection
                            print >> ofd, '      evento (t-T0) : ', TIMEnew[k], ' energia (MeV): ', Enecol[k], ' TT: ', TIMEcol[k], ' flag: ', EVcol[k]
                    if (TIMEnew[k] > t1b) and (TIMEnew[k] <t1s) or (TIMEnew[k] > t2s) and (TIMEnew[k] <t2b):
                        bkg = bkg + 1
			print >> ofd, '      fondo  (t-T0) : ', TIMEnew[k], ' energia: ', Enecol[k], ' TT: ', TIMEcol[k], ' flag: ', EVcol[k]

    print >> ofd, ''
    print >> ofd, "Source :", source
    print >> ofd, "Bkg :", bkg
    print >> ofd, ''
    print >> ofd, "N_on :", source + bkg
    print >> ofd, "N_off :", bkg
    print >> ofd, ''

    # calcolo delle significativita con il metodo di Li&Ma
    summ_src = 0
    ntot_src = 0
    summ_bkg = 0
    ntot_bkg = 0

    for i in range(len(dec_ea)):
        Vx_ea[i] = cos(dec_ea[i]) * cos(ra_ea[i])
        Vy_ea[i] = cos(dec_ea[i]) * sin(ra_ea[i])
        Vz_ea[i] = sin(dec_ea[i])

        argx_ea[i] = Vxgrb * Vx_ea[i]
        argy_ea[i] = Vygrb * Vy_ea[i]
        argz_ea[i] = Vzgrb * Vz_ea[i]
        arg_ea[i] = argx_ea[i] + argy_ea[i] + argz_ea[i]

        DELTA_ea[i] = rad2deg(arccos(arg_ea[i])) # angular distance between Earth and Source
        if (DELTA_ea[i] - ea_th + raggio) / (2. * raggio) < 0.:
            expo[i] = 0.   # ring witihin Earth region
        elif (DELTA_ea[i] - ea_th + raggio) / (2. * raggio) > 1.:
            expo[i] = 1.   # ring outside Earth region
        else:
            expo[i] = (DELTA_ea[i] - ea_th + raggio) / (2. * raggio) # VF is this true?

        if (fov - OFF[i] + raggio) / (2. * raggio) < 0.:
            expo2[i] = 0. # the ring is outside the fov
        elif (fov - OFF[i] + raggio) / (2. * raggio) > 1.:
            expo2[i] = 1. # the ring is within the fov
        else:
            expo2[i] = (fov - OFF[i] + raggio) / (2. * raggio) # is it true?

        if (((ra_punt[i] < 0) == 0) * ((ra_punt[i] > 0) == 0)): # ??????????????
            expo2[i] = 0.
        if (livetime[i] != 0) and (phase_log[i] != 1):
            expo1[i] = livetime[i] / 100. # ??????????????????
        else:
            expo1[i] = 0.

        exposure[i] = expo[i] * expo1[i] * expo2[i]

        if (TIMEnew_log[i] > t1s) and (TIMEnew_log[i] < t2s):
            summ_src = exposure[i] + summ_src
            ntot_src = ntot_src + 1

        if (TIMEnew_log[i] > t1b) and (TIMEnew_log[i] <t1s) or (TIMEnew_log[i] > t2s) and (TIMEnew_log[i] < t2b):
            summ_bkg = exposure[i] + summ_bkg
            ntot_bkg = ntot_bkg+1

    print >> ofd, "mean src not occulted ", summ_src / ntot_src
    print >> ofd, "mean bkg not occulted", summ_bkg / ntot_bkg
    mean_src = summ_src / ntot_src
    mean_bkg = summ_bkg / ntot_bkg

    if (t1s >= t2b) or (t2s <= t1b):
        tback = t2b - t1b
    else:
        if (t1s >= t1b) and (t2s <= t2b):
            tback = t2b - t1b - (t2s - t1s)
        else:
            print >> ofd, "   !!!! Scegli un altro intervallo di background !!!!"
            tback = 0

    alp = ((t2s-t1s) * mean_src) / (tback * mean_bkg)
    alp1 = alp / (1 + alp)
    alp2 = alp + 1
    print >> ofd, "source", source
    print >> ofd, "bkg", bkg / ((tback) * mean_bkg) * ((t2s - t1s) * mean_src), "(", bkg, ")"
    source1 = float(source)
    bkg1 = float(bkg)

    if ((source > 0) and (bkg > 0)):
        L1 = pow(((source1 + bkg1) / source1) * alp1, source1)
        L2 = pow(((bkg1 + source1) / bkg1) / alp2, bkg1)
        L = L1 * L2
        #print "L", alp2
        print >> ofd, 'Alpha: ', alp

        S = sqrt(-2. * log(L))
        print >> ofd, "Li&Ma sigma", S
        print >> sigmafd, S
    else:
        print >> ofd, 'Alpha: 0'
        print >> ofd, "Li&Ma sigma 0"
        print >> sigmafd, "0"

    hdulist.close()
    if ((source > 0) and (bkg > 0)):
        return alp, S, t1s, t2s, t1b, t2b
    else:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

if __name__ == '__main__':
    if (len(argv) < 11):
        print ("Wrong number of parameters.")
        print ("Usage: " + argv[0] + " GRB_l, GRB_b, radius, time_tt, t1s, t2s, t1b, t2b, fov, albedo")
        exit(0)

    l = float(argv[1])
    b = float(argv[2])
    radius = float(argv[3])
    time_tt = float(argv[4])
    t1s = float(argv[5])
    t2s = float(argv[6])
    t1b = float(argv[7])
    t2b = float(argv[8])
    fov = float(argv[9])
    albedo = float(argv[10])
    print ("Run with parameters: " + str(l) + " " + str(b) + " " + str(radius) + " " + str(time_tt)
           + " " + str(t1s) + " " + str(t2s) + " " + str(t1b) + " " + str(t2b) + " " + str(radius)
           + " " + str(fov) + " " + str(albedo))

    t1 = min(time_tt, time_tt + t1s, time_tt + t1b)
    t2 = max(time_tt, time_tt + t2s, time_tt + t2b)
    tf1 = tempfile.NamedTemporaryFile()
    log_file = tf1.name
    select_log_cmd = "$AGILE/bin/AG_select5 " + log_file + " " + log_index + " None " + str(t1) + " " + str(t2)
    print ("select log command: " + select_log_cmd)
    os.system(select_log_cmd)

    tf2 = tempfile.NamedTemporaryFile()
    evt_file = tf2.name
    select_evt_cmd = "$AGILE/bin/AG_select5 " + evt_file + " " + evt_index + " None " + str(t1) + " " + str(t2)
    print ("select evt command: " + select_evt_cmd)
    os.system(select_evt_cmd)

    out_name = "grb_analysis.out"
    sigma_name = "grb_analysis.sigma"

    ofd = open(out_name, 'w')
    sigmafd = open(sigma_name, 'w')
    ra, dec = fromGaltoRaDec(l, b)
    grb_pipe(ofd, sigmafd, evt_file, log_file, time_tt, l, b, ra, dec, t1s, t2s, t1b, t2b, radius, fov, albedo)
    ofd.close()
    sigmafd.close()
    print ("Result files: " + out_name + ", " + sigma_name)
