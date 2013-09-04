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
    if not os.path.exists('tracks/all_f'):
        os.makedirs('tracks/all_f')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/all_f'):
        os.makedirs('figures/all_f')

    grid = tracpy.inout.readgrid(loc)

    # startdates = np.array([datetime(2006, 7, 11, 0, 1)])
    startdates = np.array([datetime(2006, 2, 1, 0, 1), datetime(2006, 7, 1, 0, 1)])
    # pdb.set_trace()

    # loop through state dates
    for startdate in startdates:

        date = startdate

        # initialize counter for number of hours to increment through simulation by
        nh = 0

        # Start from the beginning and add days on for loop
        # keep running until we hit the next month
        while date.month < startdate.month+1:

            name = 'all_f/' + date.isoformat()[0:13] 

            # # If the particle trajectories have not been run, run them
            # if not os.path.exists('tracks/' + name + '.nc'):

            #     # Read in simulation initialization
            #     nstep, ndays, ff, tseas, ah, av, lon0, lat0, z0, zpar, do3d, doturb, \
            #             grid, dostream, N, T0, U, V = init.all_f(date, loc, grid=grid)

            #     # Run tracpy
            #     lonp, latp, zp, t, grid, T0, U, V \
            #         = tracpy.run.run(loc, nstep, ndays, ff, date, tseas, ah, av, \
            #                             lon0, lat0, z0, zpar, do3d, doturb, name, \
            #                             grid=grid, dostream=dostream, T0=T0, U=U, V=V)

            # # If basic figures don't exist, make them
            # if not os.path.exists('figures/' + name + '*.png'):

            # Read in and plot tracks
            d = netCDF.Dataset('tracks/' + name + '.nc')
            lonp = d.variables['lonp'][:]
            latp = d.variables['latp'][:]
            # tracpy.plotting.tracks(lonp, latp, name, grid=grid)
            # tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
            d.close()
            # Do transport plot
            tracpy.plotting.transport(name='all_f', fmod=date.isoformat()[0:13], 
                extraname=date.isoformat()[0:13], 
                Title='Transport on Shelf, for a week from ' + date.isoformat()[0:13], dmax=1.0)

            # Increment by 24 hours for next loop, to move through more quickly
            nh = nh + 24
            date = startdate + timedelta(hours=nh)
   
        # Do transport plot
        tracpy.plotting.transport(name='all_f', fmod=startdate.isoformat()[0:7], 
            extraname=startdate.isoformat()[0:7], Title='Transport on Shelf', dmax=1.0)


if __name__ == "__main__":
    run()    