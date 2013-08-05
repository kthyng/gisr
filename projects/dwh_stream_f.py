'''
Script to run the forward Deepwater Horizon drifter simulation. Start the 
drifters near Deepwater Horizon and run them forward.
This is particular to running simulations for calculating Lagrangian streamfunctions.
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
from mpl_toolkits.basemap import Basemap
import op

def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/dwh_stream_f'):
        os.makedirs('tracks/dwh_stream_f')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/dwh_stream_f'):
        os.makedirs('figures/dwh_stream_f')

    # Parameters to be rotated through
    years = np.array([2010])
    # oil spill from April 20 - July 15, 2010
    startdate = datetime(years[0], 4, 20, 0, 1)

    # Number of drifters to use
    Ns = np.array([100])#, 1000, 10000])
    # Ns = np.array([10, 20, 50, 100, 1000, 10000])

    rundays = 87 # to cover oil spill time period

    # # Do one initialization here to save grid
    # _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, grid = init.dwh_f()

    # # loop through numbers of drifters for test
    # for N in Ns:

    #     # loop through start dates for drifters
    #     for n in xrange(rundays):

    #         # Run for each output during the day too
    #         for nh in range(0,24,4):

    #             # Date for this loop
    #             date = startdate + timedelta(days=n) + timedelta(hours=nh)

    #             # Read in simulation initialization
    #             loc, nstep, ndays, ff, date, tseas, ah, av, \
    #                     lon0, lat0, z0, zpar, do3d, doturb, \
    #                     name, grid, dostream, T0, \
    #                     U, V = init.dwh_stream_f(date, N, grid=grid)

    #             # If the particle trajectories have not been run, run them
    #             if not os.path.exists('tracks/' + name + '.nc'):
    #                 lonp, latp, zp, t, grid, T0, U, V = tracpy.run.run(loc, nstep, ndays, \
    #                                                 ff, date, tseas, ah, av, \
    #                                                 lon0, lat0, z0, zpar, do3d, \
    #                                                 doturb, name, grid=grid, \
    #                                                 dostream=dostream, T0=T0, \
    #                                                 U=U, V=V)
    #             elif not os.path.exists('figures/' + name + 'tracks.png') or \
    #                  not os.path.exists('figures/' + name + 'histhexbin.png'):
    #                 d = netCDF.Dataset('tracks/' + name + '.nc')
    #                 lonp = d.variables['lonp'][:]
    #                 latp = d.variables['latp'][:]
    #                 T0 = d.variables['T0'][:]
    #                 U = d.variables['U'][:]
    #                 V = d.variables['V'][:]
    #                 tracpy.plotting.tracks(lonp, latp, name, grid=grid)
    #                 tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
    
    # Do transport plot
    tracpy.plotting.transport(name='dwh_stream_f', Title='Deepwater Horizon Spill Transport',
        fmod='*N100', dmax=2.0, N=8, llcrnrlon=-93.5, llcrnrlat=27.2, urcrnrlat=30.7,
        colormap='gray_r')

if __name__ == "__main__":
    run()