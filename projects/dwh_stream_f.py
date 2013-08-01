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
    Ns = np.array([100, 1000, 10000])
    # Ns = np.array([10, 20, 50, 100, 1000, 10000])

    rundays = 87 # to cover oil spill time period

    # Do one initialization here to save grid
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, grid = init.dwh_f()

    # loop through numbers of drifters for test
    for N in Ns:

        # loop through start dates for drifters
        for n in xrange(rundays):

            # Run for each output during the day too
            for nh in range(0,24,4):

                # Date for this loop
                date = startdate + timedelta(days=n) + timedelta(hours=nh)

                # Read in simulation initialization
                loc, nstep, ndays, ff, date, tseas, ah, av, \
                        lon0, lat0, z0, zpar, do3d, doturb, \
                        name, grid, dostream, T0, \
                        Urho, Vrho = init.dwh_stream_f(date, N, grid=grid)

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

    # ## Plot Lagrangian stream functions
    # # Which files to read in
    # Files = glob.glob('tracks/dwh_stream_f/*.nc')
    # # Files = glob.glob('tracks/bara_stream_b/*.nc')
    # Files.sort()

    # for i, File in enumerate(Files):
    #     d = netCDF.Dataset(File)
    #     if i == 0: # initialize U and V transports from first file
    #         Urho = d.variables['Urho'][:]
    #         Vrho = d.variables['Vrho'][:]
    #     else: # add in transports from subsequent simulations
    #         Urho = Urho + d.variables['Urho'][:]
    #         Vrho = Vrho + d.variables['Vrho'][:]

    # # Calculate lagrangian barotropic stream function
    # psi_i = np.cumsum(Vrho, axis=1)
    # psi_j = np.cumsum(Urho, axis=0)
    # psi = psi_j - psi_i

    # # Smaller basemap parameters.
    # llcrnrlon=-93.5; llcrnrlat=27.2; urcrnrlat=30.7
    # loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    # grid = tracpy.inout.readgrid(loc, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, 
    #                                 urcrnrlat=urcrnrlat)


    # # Make plot
    # plt.figure(figsize=(16.0375,   9.9125))
    # tracpy.plotting.background(grid=grid)
    # # contourf(grid['xr'], grid['yr'], psi, cmap='RdBu_r',
    # #         levels=np.linspace(-600,600,10), 
    # #         extend='both')
    # contourf(grid['xr'], grid['yr'], psi, cmap='RdBu_r',
    #         levels=np.linspace(-120000,120000,12), 
    #         extend='both')
    # colorbar()
    # # Add initial drifter location (all drifters start at the same location)
    # lon0 = d.variables['lonp'][0,0]
    # lat0 = d.variables['latp'][0,0]
    # x0, y0 = grid['basemap'](lon0, lat0)
    # plot(x0, y0, 'go', markersize=10)
    # plt.savefig('figures/dwh_stream_f/stream',bbox_inches='tight')

    # # # Compare with tracks plot and histogram
    # # d = netCDF.MFDataset('tracks/bara_stream_b/*.nc',aggdim='ntrac')
    # # lonp = d.variables['lonp'][:]
    # # latp = d.variables['latp'][:]
    # # tracpy.plotting.tracks(lonp, latp, 'bara_stream_b/overall', grid=grid)
    # # tracpy.plotting.hist(lonp, latp, 'bara_stream_b/overall', grid=grid, which='hexbin')
    # # tracpy.plotting.hist(lonp, latp, 'bara_stream_b/overall-alllocations', grid=grid, which='hexbin', tind='all')
    # # Compare with tracks plot and histogram
    # d = netCDF.MFDataset('tracks/dwh_stream_f/*.nc',aggdim='ntrac')
    # lonp = d.variables['lonp'][:]
    # latp = d.variables['latp'][:]
    # tracpy.plotting.hist(lonp, latp, 'dwh_stream_f/overall-alllocations', grid=grid, which='hexbin', tind='all')
    # tracpy.plotting.tracks(lonp, latp, 'dwh_stream_f/overall', grid=grid)
    # tracpy.plotting.hist(lonp, latp, 'dwh_stream_f/overall', grid=grid, which='hexbin')


if __name__ == "__main__":
    run()