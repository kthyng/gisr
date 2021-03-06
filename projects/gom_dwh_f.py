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

    units = 'seconds since 1858-11-18' #????
    loc = 'http://omgsrv1.meas.ncsu.edu:8080/thredds/dodsC/fmrc/sabgom/SABGOM_Forecast_Model_Run_Collection_best.ncd'

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/gom_dwh_f'):
        os.makedirs('tracks/gom_dwh_f')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/gom_dwh_f'):
        os.makedirs('figures/gom_dwh_f')

    # Parameters to be rotated through
    years = np.array([2013])
    # forecast model data available
    startdate = datetime(years[0], 7, 31, 15, 1)

    # Number of drifters to use
    Ns = np.array([100])#, 1000, 10000])
    # Ns = np.array([10, 20, 50, 100, 1000, 10000])

    rundays = 4

    # Do one initialization here to save grid
    grid = tracpy.inout.readgrid(loc)

    # loop through numbers of drifters for test
    for N in Ns:

        # loop through start dates for drifters
        for n in xrange(rundays):

            # Run for each output during the day too
            for nh in range(0,24,3):

                # Date for this loop
                date = startdate + timedelta(days=n) + timedelta(hours=nh)

                # Read in simulation initialization
                loc, nstep, ndays, ff, date, tseas, ah, av, \
                        lon0, lat0, z0, zpar, do3d, doturb, \
                        name, grid, dostream, T0, \
                        U, V = init.gom_dwh_f(date, N, grid=grid)

                # If the particle trajectories have not been run, run them
                if not os.path.exists('tracks/' + name + '.nc'):
                    lonp, latp, zp, t, grid, T0, U, V = tracpy.run.run(loc, nstep, ndays, \
                                                    ff, date, tseas, ah, av, \
                                                    lon0, lat0, z0, zpar, do3d, \
                                                    doturb, name, grid=grid, \
                                                    dostream=dostream, T0=T0, \
                                                    U=U, V=V)

                else: # if the files already exist, just read them in for plotting
                    d = netCDF.Dataset('tracks/' + name + '.nc')
                    lonp = d.variables['lonp'][:]
                    latp = d.variables['latp'][:]
                    T0 = d.variables['T0'][:]
                    U = d.variables['U'][:]
                    V = d.variables['V'][:]

                # If the particle trajectories have not been plotted, plot them
                if not os.path.exists('figures/' + name + 'tracks.png'):
                    tracpy.plotting.tracks(lonp, latp, name, grid=grid)
                if not os.path.exists('figures/' + name + 'histhexbin.png'):
                    tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')

                # # Plot tracks and histograms for these drifters
                # # but only those outside the shelf
                # fh = grid['trirllrho'].nn_interpolator(grid['h'].flatten())
                # hp = fh(lonp[:,0],latp[:,0]) #depths at starting lon/lat
                # ind = hp > 500. # want to know which drifters start outside shelf break
                # name1 = name + 'tracks_outershelf'
                # if not os.path.exists('figures/' + name1 + '.png'):
                #     tracpy.plotting.tracks(lonp[ind,:], latp[ind,:], name1, grid=grid)
                # name2 = name + 'histhexbin_outershelf'
                # if not os.path.exists('figures/' + name2 + '.png'):
                #     tracpy.plotting.hist(lonp[ind,:], latp[ind,:], name2, grid=grid, \
                #                         which='hexbin')

if __name__ == "__main__":
    run()