'''
Script to run the forward Galveston drifter simulation. Start the 
drifters near Galveston and run them forward.
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
    # Need 1 min to get correct model output out
    startdates = np.array([datetime(years[0], 5, 23, 0, 1)])#,
                            # datetime(years[0], 6, 21, 0, 1)])

    # how many days to start drifters from, starting from startdate
    rundays = 6

    grid = tracpy.inout.readgrid(loc)

    # loop through start dates for drifters
    for startdate in startdates:
        for n in xrange(rundays):
            # Run for each hour during the day too
            for nh in range(0,24,4):

                # Date for this loop
                date = startdate + timedelta(days=n) + timedelta(hours=nh)

                # Read in simulation initialization
                loc, nstep, ndays, ff, date, tseas, ah, av, \
                        lon0, lat0, z0, zpar, do3d, doturb, \
                        name, grid, dostream, T0, \
                        U, V = init.galv_f(date, grid=grid)

                # If the particle trajectories have not been run, run them
                if not os.path.exists('tracks/' + name + '.nc'):
                    lonp, latp, zp, t, grid, T0, U, V = tracpy.run.run(loc, nstep, ndays, \
                                                    ff, date, tseas, ah, av, \
                                                    lon0, lat0, z0, zpar, do3d, \
                                                    doturb, name, grid=grid, \
                                                    dostream=dostream, T0=T0, \
                                                    U=U, V=V)

                elif not os.path.exists('figures/' + name + 'tracks.png') or \
                     not os.path.exists('figures/' + name + 'histhexbin.png'):
                    d = netCDF.Dataset('tracks/' + name + '.nc')
                    lonp = d.variables['lonp'][:]
                    latp = d.variables['latp'][:]
                    T0 = d.variables['T0'][:]
                    U = d.variables['U'][:]
                    V = d.variables['V'][:]
                    tracpy.plotting.tracks(lonp, latp, name, grid=grid)
                    tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
    
    # Do transport plot
    tracpy.plotting.transport(name='galv_f', Title='Transport from Galveston', dmax=1.5)

# This is so the script can be run using reference 
# projects/[projectname].py
if __name__ == "__main__":
    run()    