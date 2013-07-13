'''
Script to run the backward Barataria drifter simulation. Start the 
drifters near Barataria and run them backward.
'''

import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
import glob

def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/bara_b'):
        os.makedirs('tracks/bara_b')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/bara_b'):
        os.makedirs('figures/bara_b')

    # Parameters to be rotated through
    years = np.array([2010])
    # oil spill from April 20 - July 15, 2010
    startdate = datetime(years[0], 4, 20, 0)

    ndays = 87 # to cover oil spill time period

    # Do one initialization here to save grid
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, grid = init.bara_b()

    # loop through start dates for drifters
    for n in xrange(ndays):

        # Date for this loop
        date = startdate + timedelta(days=n)

        # Read in simulation initialization
        loc, nstep, ndays, ff, date, tseas, ah, av, lon0, lat0, z0, \
                zpar, do3d, doturb, name, grid = init.bara_b(date, grid=grid)

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc'):
            lonp, latp, zp, t, grid = tracpy.run.run(loc, nstep, ndays, \
                                            ff, date, tseas, ah, av, \
                                            lon0, lat0, z0, zpar, do3d, \
                                            doturb, name)

        else: # if the files already exist, just read them in for plotting
            d = netCDF.Dataset('tracks/' + name + '.nc')
            lonp = d.variables['lonp'][:]
            latp = d.variables['latp'][:]


        # If the particle trajectories have not been plotted, plot them
        if not os.path.exists('figures/' + name + 'tracks.png'):
            tracpy.plotting.tracks(lonp, latp, name, grid=grid)
        if not os.path.exists('figures/' + name + 'histhexbin.png'):
            tracpy.plotting.hist(lonp, latp, name, grid=grid, \
                                which='hexbin', bins=(40,40))
