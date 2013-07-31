'''
Script to run the backward Barataria Bay drifter simulation. Start the 
drifters near Barataria Bay and run them backward.
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

def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/bara_stream_b'):
        os.makedirs('tracks/bara_stream_b')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/bara_stream_b'):
        os.makedirs('figures/bara_stream_b')

    # Parameters to be rotated through
    years = np.array([2010])
    # oil spill from April 20 - July 15, 2010
    startdate = datetime(years[0], 4, 20, 0)

    # Number of drifters to use
    Ns = np.array([100])#, 1000, 10000])
    # Ns = np.array([10, 20, 50, 100, 1000, 10000])

    rundays = 87 # to cover oil spill time period

    # Do one initialization here to save grid
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, grid = init.dwh_f()

    # loop through start dates for drifters
    for n in xrange(rundays):

        # loop through numbers of drifters for test
        for N in Ns:

            # Run for each output during the day too
            for nh in range(0,24,4):

                # Date for this loop
                date = startdate + timedelta(days=n) + timedelta(hours=nh)

                # Read in simulation initialization
                loc, nstep, ndays, ff, date, tseas, ah, av, \
                        lon0, lat0, z0, zpar, do3d, doturb, \
                        name, grid, dostream, T0, \
                        Urho, Vrho = init.bara_stream_b(date, N, grid=grid)

                # If the particle trajectories have not been run, run them
                if not os.path.exists('tracks/' + name + '.nc'):
                    lonp, latp, zp, t, grid, Urho, Vrho = tracpy.run.run(loc, nstep, ndays, \
                                                    ff, date, tseas, ah, av, \
                                                    lon0, lat0, z0, zpar, do3d, \
                                                    doturb, name, grid=grid, \
                                                    dostream=dostream, T0=T0, \
                                                    Urho=Urho, Vrho=Vrho)

                else: # if the files already exist, just read them in for plotting
                    d = netCDF.Dataset('tracks/' + name + '.nc')
                    lonp = d.variables['lonp'][:]
                    latp = d.variables['latp'][:]
                    Urho = d.variables['Urho'][:]
                    Vrho = d.variables['Vrho'][:]

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

                # # Plot Lagrangian stream functions
                # Lx = np.zeros(Urho.shape)
                # for i in xrange(1,Urho.shape[0]):
                #     Lx[i,:] = Lx[i-1,:] + Urho[i,:]
                # Ly = np.zeros(Urho.shape)
                # for j in xrange(1,Urho.shape[1]):
                #     Ly[:,j] = Ly[:,j-1] - Vrho[:,j]
                # pdb.set_trace()

if __name__ == "__main__":
    run()