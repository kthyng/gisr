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

# on pong
loc = ['/pong/raid/kthyng/forecast/roms_his_20130101_a_n0.nc','http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc']

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

def init(grid, date, lon0, lat0):

    nsteps = 25 # limiting steps 
    N = 5 # 5 outputs between model outputs
    ndays = 30 # days
    ff = 1 # forward for now but ALSO TRY BACKWARD FROM FOUND LOCATIONS. SHOULD DO BOTH!

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 5.
    av = 0. # m^2/s

    # Number of drift cards set out
    ndrift = 1770

    # startdate, startlon, startlat, enddate, endlon, endlat, idnum = reduce(grid, ndrift)

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 2

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    # Transport would be nice for pictures, and starting in a single location lends itself
    # to this nicely.
    dostream = 1

    # convert date to number
    datenum = netCDF.date2num(date, units)

    # Number of model outputs to use
    tout = np.int((ndays*(24*3600))/tseas)

    # Figure out what files will be used for this tracking - to get tinds for
    # the following calculation
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)

    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()

    # DO THIS ONCE I KNOW THE POSITION
    Ueast = uf[1:, 1:-1, 0].copy()
    ind = Ueast<0 # positive here is outflowing
    Ueast[ind] = 0 # don't want to count negative values of Ueast
    Uwest = uf[:-1, 1:-1, 0].copy()
    ind = Uwest>0 # negative is outflowing 
    Uwest[ind] = 0
    Vnorth = vf[1:-1, 1:, 0].copy()
    ind = Vnorth<0 # positive is outflowing
    Vnorth[ind] = 0
    Vsouth = vf[1:-1, :-1, 0].copy()
    ind = Vsouth>0 # negative is outflowing from cell
    Vsouth[ind] = 0
    Tgrid = abs(Ueast) + abs(Vnorth) + abs(Uwest) + abs(Vsouth)

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

    grid = tracpy.inout.readgrid(loc)

    # Read in deployment info
    ndrift = 3490 # number of driftcards deployed
    startdate, startlon, startlat, enddate, endlon, endlat, idnum = read(ndrift)

    # Eliminate those outside the numerical domain
    lon0, lat0 = tracpy.tools.check_points(lon0, lat0, grid)

    # Loop through deployments

        # Check for in domain and if not, continue

        # Initialize parameters
        startdate, startlon, startlat, enddate, endlon, endlat, idnum, \
                nsteps, ndays, ff, tseas, ah, av, ndrift, z0, zpar, \
                do3d, doturb, dostream = init(grid)

        # start drifters in chunks of exact same date/time
        # idrift = 0
        # tdrift = startdate[idrift] #drifter time
        idrifts = find(startdate==startdate[find(~np.isnan(startdate))[0]]) # all drifters with same start date/time
        lon0 = startlon[idrifts]
        lat0 = startlat[idrifts]
        date = netCDF.num2date(startdate[idrifts][0], units)

        pdb.set_trace()

        name = 'driftcards/' + date.isoformat()

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