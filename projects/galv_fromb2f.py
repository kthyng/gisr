'''
Script to run from the backward Galveston drifter simulation origins, forward.
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

units = 'seconds since 1970-01-01'

def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/galv_fromb2f'):
        os.makedirs('tracks/galv_fromb2f')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/galv_fromb2f'):
        os.makedirs('figures/galv_fromb2f')

    # Parameters to be rotated through
    years = np.array([2010])
    # oil spill from April 20 - July 15, 2010
    # startdate = datetime(years[0], 4, 20, 0, 1)
    # Interesting time periods where backward drifters cross the shelf:
    # 5-23-10 to 5-28-10 and 6-21-10 to 6-26-10
    # Need 1 min to get correct model output out
    startdates = np.array([datetime(years[0], 5, 23, 0, 1)])#,
                            # datetime(years[0], 6, 21, 0, 1)])

    # how many days to start drifters from, starting from startdate
    rundays = 6

    # Do one initialization here to save grid
    _, _, _, _, _, _, _, _, _, _, _, _, grid = init.galv_fromb2f()

    # loop through start dates for drifters
    # pdb.set_trace()
    Files = glob.glob('tracks/galv_b/2010-05-*')
    for File in Files:

        d = netCDF.Dataset(File)
        tp = d.variables['tp'][:]
        backward_startdate = netCDF.num2date(tp[0],units)

        # calculate last time step (since wasn't saved correctly), which is start date
        date = backward_startdate-timedelta(minutes=48)*(len(tp)-1)

        # get final positions, which are start positions for this run
        lon0 = d.variables['lonp'][:,-1]
        lat0 = d.variables['latp'][:,-1]
        # don't check points because they are already within the domain from
        # running backward
        # lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)

        # Read in simulation initialization
        loc, nstep, ndays, ff, tseas, ah, av, z0, \
                zpar, do3d, doturb, name, grid = init.galv_fromb2f(date, grid=grid)

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
            tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')

            #     # Save up tracks to plot together by month
            #     if n == 0:
            #         lonpsavem = lonp
            #         latpsavem = latp
            #     else:
            #         lonpsavem = np.vstack((lonpsavem,lonp))
            #         latpsavem = np.vstack((latpsavem,latp))
            #     if n == 0: # initially, give month a value
            #         month = date.month
            #     # want to know when we are at a new month or the end of the year
            #     elif (date.month != month) or \
            #             (date.month == 12 and date.day == 31): 
            #         # First plot previous month
            #         # take off day and month and add on previous month instead
            #         name = name[:-5] + str(month).zfill(2)
            #         tracpy.plotting.tracks(lonpsavem, latpsavem, name, grid=grid)
            #         tracpy.plotting.hist(lonpsavem, latpsavem, name, grid=grid, \
            #                             which='hexbin', bins=(40,40))
            #         # Reset month to next month value
            #         month = date.month
            #         # Save month arrays for year plot
            #         if month == 1:
            #             lonpsavey = lonp
            #             latpsavey = latp
            #         else:
            #             lonpsavey = np.vstack((lonpsavey,lonpsavem))
            #             latpsavey = np.vstack((latpsavey,latpsavem))
            #         # Reset save arrays for month plots
            #         lonpsavem = lonp
            #         latpsavem = latp

            # # Plot year summaries
            # name = str(year)
            # tracpy.plotting.tracks(lonpsavey, latpsavey, name, grid=grid)
            # tracpy.plotting.hist(lonpsavey, latpsavey, name, grid=grid, \
            #                     which='hexbin', bins=(40,40))

            # ## Weatherband plotting
            # # Read in all tracks
            # files = np.sort(glob.glob('tracks/galv_b/*.nc')) # sorted list of file names

    # # Plot forward and backward together for Galveston
    # # May time period: 5/23/10-5/28/10
    # db = netCDF.MFDataset('tracks/galv_b/2010-05-*.nc',aggdim='ntrac')
    # lonb = np.fliplr(db.variables['lonp'][:]) # flip to be forward in time
    # latb = np.fliplr(db.variables['latp'][:]) # flip to be forward in time
    # df = netCDF.MFDataset('tracks/galv_b/2010-06-*.nc',aggdim='ntrac')
    # lonf = df.variables['lonp'][:]
    # latf = df.variables['latp'][:]
    # # Combined
    # lonp = np.concatenate((lonb,lonf),axis=1)
    # latp = np.concatenate((latb,latf),axis=1)
    # # Plot
    # name = '2010-05-23-back+forward'
    # tracpy.plotting.tracks(lonp, latp, name, grid=grid)
    # tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
                                            
    # June time period: 6/21/10-6/26/10


if __name__ == "__main__":
    run()    