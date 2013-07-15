'''
Script to run the forward Galveston drifter simulation. Start the 
drifters near Galveston and run them forward.
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
    if not os.path.exists('tracks/galv_f'):
        os.makedirs('tracks/galv_f')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/galv_f'):
        os.makedirs('figures/galv_f')

    # Parameters to be rotated through
    years = np.array([2010])
    # Interesting time periods where backward drifters cross the shelf:
    # 5-23-10 to 5-28-10 and 6-21-10 to 6-26-10
    startdates = np.array([datetime(years[0], 5, 23, 0),
                            datetime(years[0], 6, 21, 0)])

    # how many days to start drifters from, starting from startdate
    rundays = 6

    # Do one initialization here to save grid
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, grid = init.galv_f()

    # loop through start dates for drifters
    for startdate in startdates:
        for n in xrange(rundays):

            # Date for this loop
            date = startdate + timedelta(days=n)

            # Read in simulation initialization
            loc, nstep, ndays, ff, date, tseas, ah, av, lon0, lat0, z0, \
                    zpar, do3d, doturb, name, grid = init.galv_f(date, grid=grid)

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
                                    which='hexbin')

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

# This is so the script can be run using reference 
# projects/[projectname].py
if __name__ == "__main__":
    run()    