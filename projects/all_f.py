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

npieces = 1# 20 # number of pieces to divide starting locations for drifters into, in x direction

def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/all_f/'):
        os.makedirs('tracks/all_f/')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/all_f/'):
        os.makedirs('figures/all_f/')
        
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    grid = tracpy.inout.readgrid(loc)

    overallstartdate = datetime(2004, 1, 1, 0, 1)
    overallstopdate = datetime(2005, 1, 1, 0, 1)

    # startdates = np.array([datetime(2004, 2, 1, 0, 1), datetime(2004, 7, 1, 0, 1),
    #                         datetime(2005, 2, 1, 0, 1), datetime(2005, 7, 1, 0, 1),
    #                         datetime(2006, 2, 1, 0, 1), datetime(2006, 7, 1, 0, 1),
    #                         datetime(2007, 2, 1, 0, 1), datetime(2007, 7, 1, 0, 1),
    #                         datetime(2008, 2, 1, 0, 1), datetime(2008, 7, 1, 0, 1),
    #                         datetime(2009, 2, 1, 0, 1), datetime(2009, 7, 1, 0, 1),
    #                         datetime(2010, 2, 1, 0, 1), datetime(2010, 7, 1, 0, 1),
    #                         datetime(2011, 2, 1, 0, 1), datetime(2011, 7, 1, 0, 1)])

    # # loop through state dates
    # for startdate in startdates:

    date = overallstartdate

    # initialize counter for number of hours to increment through simulation by
    # nh = 0

    # Start from the beginning and add days on for loop
    # keep running until we hit the next month
    while date < overallstopdate:

        # Loop over the domain, divided into pieces to make files more manageable
        for i in xrange(npieces):

            name = 'all_f/' + date.isoformat()[0:13] + '_' + str(i)

            # If the particle trajectories have not been run, run them
            if not os.path.exists('tracks/' + name + '.nc'):

                # Read in simulation initialization
                nstep, N, ndays, ff, tseas, ah, av, lon0, lat0, z0, zpar, do3d, doturb, \
                        grid, dostream, T0, U, V = init.all_f(date, loc, 
                                                    npiece=[i/float(npieces), (i+1)/float(npieces)], grid=grid)
                # pdb.set_trace()
                # Run tracpy
                lonp, latp, zp, t, grid, T0, U, V \
                    = tracpy.run.run(loc, nstep, ndays, ff, date, tseas, ah, av, \
                                        lon0, lat0, z0, zpar, do3d, doturb, name, N=N,  \
                                        grid=grid, dostream=dostream, T0=T0, U=U, V=V)

            # # If basic figures don't exist, make them
            # if not os.path.exists('figures/' + name + '*.png'):

                # Read in and plot tracks
                d = netCDF.Dataset('tracks/' + name + '.nc')
                lonp = d.variables['lonp'][:]
                latp = d.variables['latp'][:]
                # tracpy.plotting.tracks(lonp, latp, name, grid=grid)
                # tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
                d.close()
                # # Do transport plot
                # tracpy.plotting.transport(name='all_f/N=5_dx=8/25days', fmod=date.isoformat()[0:13], 
                #     extraname=date.isoformat()[0:13], 
                #     Title='Transport on Shelf, for a week from ' + date.isoformat()[0:13], dmax=1.0)

            # Increment by 24 hours for next loop, to move through more quickly
            # nh = nh + 24
        date = date + timedelta(hours=24*7)
   
        # # Do transport plot
        # tracpy.plotting.transport(name='all_f/N=5_dx=8/25days', fmod=startdate.isoformat()[0:7] + '*', 
        #     extraname=startdate.isoformat()[0:7], Title='Transport on Shelf', dmax=1.0)


if __name__ == "__main__":
    run()    
