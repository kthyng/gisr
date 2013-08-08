'''
Script to run drifters backward from Galveston Bay to examine the Bay's 
connectivity with the shelf region.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
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

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/galvcon_b'):
        os.makedirs('tracks/galvcon_b')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/galvcon_b'):
        os.makedirs('figures/galvcon_b')

    grid = tracpy.inout.readgrid(loc)

    # Start from the beginning and add days on for loop
    date = datetime(2004, 1, 1, 0, 1)

    # initialize counter for number of hours to increment through simulation by
    nh = 0

    # keep running until we hit 2012
    while date.year < 2012:

        name = 'galvcon_b/' + date.isoformat()[0:13] 

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc'):

            # Read in simulation initialization
            nstep, ndays, ff, tseas, ah, av, lon0, lat0, z0, zpar, do3d, doturb, \
                    grid, dostream, N, T0, U, V = init.galvcon_b(date, loc, grid=grid)

            # Run tracpy
            lonp, latp, zp, t, grid, T0, U, V \
                = tracpy.run.run(loc, nstep, ndays, ff, date, tseas, ah, av, \
                                    lon0, lat0, z0, zpar, do3d, doturb, name, \
                                    grid=grid, dostream=dostream, T0=T0, U=U, V=V)

        # If basic figures don't exist, make them
        if not os.path.exists('figures/' + name + '*.png'):

            # Read in and plot tracks
            d = netCDF.Dataset('tracks/' + name + '.nc')
            lonp = d.variables['lonp'][:]
            latp = d.variables['latp'][:]
            tracpy.plotting.tracks(lonp, latp, name, grid=grid)
            tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')

        # Increment by 4 hours for next loop
        nh = nh + 4
        date = date + timedelta(hours=nh)
   
    # Do transport plot
    tracpy.plotting.transport(name='galvcon_b', Title='Transport to Galveston', dmax=1.5)


if __name__ == "__main__":
    run()    