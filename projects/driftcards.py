from __future__ import division
import csv
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
from matplotlib.mlab import find, Path

units = 'seconds since 1970-01-01'

# Location of TXLA model output
loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

def read(ndrift):

    # Read in all initial drift card locations
    idnum = np.zeros(ndrift).astype(int)
    startdate = np.zeros(ndrift)*np.nan
    startlat = np.zeros(ndrift)*np.nan
    startlon = np.zeros(ndrift)*np.nan

    with open('projects/driftcards/Deployment_data.txt','rU') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        for i, row in enumerate(reader):
            idnum[i] = row[0]
            startlat[i] = row[2]
            startlon[i] = row[3]

            temp = row[1]
            year = int(temp.split()[0].split('-')[0])
            month = int(temp.split()[0].split('-')[1])
            day = int(temp.split()[0].split('-')[2])
            hour = int(temp.split()[1].split(':')[0])
            if temp.split()[2] == 'PM' and hour is not 12:
                hour = hour + 12
            minutes = int(temp.split()[1].split(':')[1])
            startdate[i] = netCDF.date2num(datetime(year, month, day, hour, minutes), units)

    # ind = idnum.argsort()
    # Find end locations for given driftcards
    enddate = np.zeros(ndrift)*np.nan
    endlat = np.zeros(ndrift)*np.nan
    endlon = np.zeros(ndrift)*np.nan

    with open('projects/driftcards/Recovery_data.txt','rU') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        for i, row in enumerate(reader):
            # this is the index that corresponds to the deployment data
            # can use this to keep everything referenced together
            ind = find(idnum==int(row[0]))
            endlat[ind] = row[2]
            endlon[ind] = row[3]

            temp = row[1]
            year = int(temp.split()[0].split('-')[0])
            month = int(temp.split()[0].split('-')[1])
            day = int(temp.split()[0].split('-')[2])
            hour = int(temp.split()[1].split(':')[0])
            if temp.split()[2] == 'PM' and hour is not 12:
                hour = hour + 12
            minutes = int(temp.split()[1].split(':')[1])
            enddate[ind] = netCDF.date2num(datetime(year, month, day, hour, minutes), units)

    # All arrays are referenced to the same ordering of the drifter idnum

    return startdate, startlon, startlat, enddate, endlon, endlat, idnum

def reduce(grid, ndrift):

    startdate, startlon, startlat, enddate, endlon, endlat, idnum = read(ndrift)

    # If covering the whole domain, need to exclude points outside domain.
    # Use info just inside domain so points aren't right at the edge.
    xvert = np.hstack((np.flipud(grid['lonr'][1,:]),grid['lonr'][:,1],
        grid['lonr'][-2,:],np.flipud(grid['lonr'][:,-2])))
    yvert = np.hstack((np.flipud(grid['latr'][1,:]),grid['latr'][:,1],
        grid['latr'][-2,:],np.flipud(grid['latr'][:,-2])))
    verts = np.vstack((xvert,yvert))
    # Form path
    path = Path(verts.T)
    # nan out all drifters that do not start within the TXLA model domain
    for i in xrange(ndrift):
        if not path.contains_point(np.vstack((startlon[i], startlat[i]))):
            startlon[i] = np.nan
            startlat[i] = np.nan
            startdate[i] = np.nan
            endlon[i] = np.nan
            endlat[i] = np.nan

    # nan out drifters that have not been found yet
    for i in xrange(ndrift):
        if np.isnan(enddate[i]):
            startlon[i] = np.nan
            startlat[i] = np.nan
            startdate[i] = np.nan
  
    return startdate, startlon, startlat, enddate, endlon, endlat, idnum

def init(grid):

    nsteps = 5 # 5 time interpolation steps
    ndays = 90#180 # in days, about 3 months
    ff = 1 # forward for now but ALSO TRY BACKWARD FROM FOUND LOCATIONS. SHOULD DO BOTH!

    # FINISH INIT TO SEE ABOUT GETTING A RUN TO WORK

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 20.
    av = 0. # m^2/s

    # Number of drift cards set out
    ndrift = 1770

    startdate, startlon, startlat, enddate, endlon, endlat, idnum = reduce(grid, ndrift)

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 0

    return startdate, startlon, startlat, enddate, endlon, endlat, idnum, \
            nsteps, ndays, ff, tseas, ah, av, ndrift, z0, zpar, do3d, doturb, dostream


def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/driftcards'):
        os.makedirs('tracks/driftcards')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/driftcards'):
        os.makedirs('figures/driftcards')
 
    # Find drifters that started within the TX-LA numerical domain
    # Use path to do so
    grid = tracpy.inout.readgrid(loc)
    # lonr = grid['lonr']
    # latr = grid['latr']

    # Initialize parameters
    startdate, startlon, startlat, enddate, endlon, endlat, idnum, \
            nsteps, ndays, ff, tseas, ah, av, ndrift, z0, zpar, \
            do3d, doturb, dostream = init(grid)

    # For drifters that started within domain and have been found, order by
    # start date and time
    # Start drifters at drift card locations, and have a different simulation for
    # each set of drifters by start date/time
    # Run for three months?

    # Don't need to go in order I suppose
    # # Sort arrays by startdate
    # startind = np.argsort(startdate) #these indices sort arrays by startdate

    # THIS NEEDS TO BE IN SOME SORT OF LOOP

    # start drifters in chunks of exact same date/time
    # idrift = 0
    # tdrift = startdate[idrift] #drifter time
    idrifts = find(startdate==startdate[find(~np.isnan(startdate))[0]]) # all drifters with same start date/time
    lon0 = startlon[idrifts]
    lat0 = startlat[idrifts]
    date = netCDF.num2date(startdate[idrifts][0], units)

    pdb.set_trace()

    name = 'driftcards/' + date.isoformat()

    # COULD RUN MORE DRIFTERS THAN THE NUMBER OF DRIFT CARDS

    # Don't have the output available for these runs yet

    # Run tracpy
    lonp, latp, zp, t, grid, T0, U, V \
        = tracpy.run.run(loc, nsteps, ndays, ff, date, tseas, ah, av, \
                            lon0, lat0, z0, zpar, do3d, doturb, name, \
                            grid=grid, dostream=dostream)


    # nan out drifters once they have been used to run tracpy
    startlon[idrifts] = np.nan
    startlat[idrifts] = np.nan
    startdate[idrifts] = np.nan

    pdb.set_trace()

if __name__ == "__main__":
    run()    