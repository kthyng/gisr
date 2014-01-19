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

def read(ndrift, fileloc):

    # Read in all initial drift card locations
    idnum = np.zeros(ndrift).astype(int)
    startdate = np.zeros(ndrift)*np.nan
    startlat = np.zeros(ndrift)*np.nan
    startlon = np.zeros(ndrift)*np.nan

    with open(fileloc,'rU') as csvfile:
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

    # # ind = idnum.argsort()
    # # Find end locations for given driftcards
    # enddate = np.zeros(ndrift)*np.nan
    # endlat = np.zeros(ndrift)*np.nan
    # endlon = np.zeros(ndrift)*np.nan

    # with open('projects/driftcards/Recovery_data.txt','rU') as csvfile:
    #     reader = csv.reader(csvfile,delimiter=',')
    #     for i, row in enumerate(reader):
    #         # this is the index that corresponds to the deployment data
    #         # can use this to keep everything referenced together
    #         ind = find(idnum==int(row[0]))
    #         endlat[ind] = row[2]
    #         endlon[ind] = row[3]

    #         temp = row[1]
    #         year = int(temp.split()[0].split('-')[0])
    #         month = int(temp.split()[0].split('-')[1])
    #         day = int(temp.split()[0].split('-')[2])
    #         hour = int(temp.split()[1].split(':')[0])
    #         if temp.split()[2] == 'PM' and hour is not 12:
    #             hour = hour + 12
    #         minutes = int(temp.split()[1].split(':')[1])
    #         enddate[ind] = netCDF.date2num(datetime(year, month, day, hour, minutes), units)

    # All arrays are referenced to the same ordering of the drifter idnum

    return startdate, startlon, startlat, idnum
    # return startdate, startlon, startlat, enddate, endlon, endlat, idnum

def readThenReduce(grid, ndrift, fileloc='projects/driftcards/Deployment_data.txt'):

    startdate, startlon, startlat, idnum = read(ndrift, fileloc)
    # startdate, startlon, startlat, enddate, endlon, endlat, idnum = read(ndrift)

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
            # endlon[i] = np.nan
            # endlat[i] = np.nan

    # # nan out drifters that have not been found yet
    # for i in xrange(ndrift):
    #     if np.isnan(enddate[i]):
    #         startlon[i] = np.nan
    #         startlat[i] = np.nan
    #         startdate[i] = np.nan
    # pdb.set_trace()
    # Eliminate nan'ed entries (only save non-nan'ed entries)
    ind = ~np.isnan(startlon)
    startlon = startlon[ind]
    startlat = startlat[ind]
    startdate = startdate[ind]
    idnum = idnum[ind]
  
    return startdate, startlon, startlat, idnum

def init(grid, datenum, lon0, lat0):

    nsteps = 25 # limiting steps 
    N = 5 # 5 outputs between model outputs
    ndays = 30 # days
    ff = 1 # forward for now but ALSO TRY BACKWARD FROM FOUND LOCATIONS. SHOULD DO BOTH!

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 5.
    av = 0. # m^2/s

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

    # # convert date to number
    # datenum = netCDF.date2num(date, units)

    # Number of model outputs to use
    tout = np.int((ndays*(24*3600))/tseas)

    # Figure out what files will be used for this tracking - to get tinds for
    # the following calculation
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)

    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0], grid, nc, z0, zpar)
    nc.close()

    # Interpolate to get starting position in grid space
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    xstart0 = np.ceil(xstart0[0]); ystart0 = ystart0[0] # they are all the same
    Ueast = uf[xstart0, ystart0, 0]
    if Ueast < 0: Ueast=0. # positive here is outflowing and don't want to count negative values of Ueast
    Uwest = uf[xstart0-1, ystart0, 0]
    if Uwest > 0: Uwest=0. # negative here is outflowing and don't want to count positive values of Ueast
    Vnorth = vf[xstart0, ystart0, 0]
    if Vnorth < 0: Vnorth=0. # positive here is outflowing and don't want to count negative values of Ueast
    Vsouth = vf[xstart0, ystart0, 0]
    if Vsouth > 0: Vsouth=0. # negative here is outflowing and don't want to count positive values of Ueast
    # drifter initial volume transport is divided by the number of drifters being run
    T0 = np.ones(lon0.size)*((abs(Ueast) + abs(Vnorth) + abs(Uwest) + abs(Vsouth))/lon0.size)

    # Initialize the arrays to save the transports on the grid in the loop.
    # These arrays aggregate volume transport when a drifter enters or exits a grid cell
    # These should start at zero since we don't know which way things will travel yet
    U = np.ma.zeros(grid['xu'].shape,order='F')
    V = np.ma.zeros(grid['xv'].shape,order='F')

    return nsteps, N, ndays, ff, tseas, ah, av, z0, zpar, do3d, doturb, dostream, T0, U, V


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

    ndriftmodel = 1000 # do 1000 drifters for each set of driftcards

    # Read in deployment info and eliminate outside domain
    ndrift = 3490 # number of driftcards deployed
    # These are drift cards that start in numerical domain
    startdate, startlon, startlat, idnum = readThenReduce(grid, ndrift)

    # Find unique start dates and indices for start dates
    dates, idates = np.unique(startdate, return_index=True)

    # Loop through deployments
    for i, date in enumerate(dates):

        idate = idates[i] # Index in arrays
        lon0 = np.ones(ndriftmodel)*startlon[idate]
        lat0 = np.ones(ndriftmodel)*startlat[idate]

        # pdb.set_trace()

        # Initialize parameters
        nsteps, N, ndays, ff, tseas, ah, av, z0, zpar, \
            do3d, doturb, dostream, T0, U, V = init(grid, date, lon0, lat0)

        # pdb.set_trace()

        date = netCDF.num2date(date, units)

        # # start drifters in chunks of exact same date/time
        # # idrift = 0
        # # tdrift = startdate[idrift] #drifter time
        # idrifts = find(startdate==startdate[find(~np.isnan(startdate))[0]]) # all drifters with same start date/time
        # lon0 = startlon[idrifts]
        # lat0 = startlat[idrifts]
        # date = netCDF.num2date(startdate[idrifts][0], units)

        # pdb.set_trace()

        name = 'driftcards/' + date.isoformat()

        # Run tracpy
        lonp, latp, zp, t, grid, T0, U, V \
            = tracpy.run.run(loc, nsteps, ndays, ff, date, tseas, ah, av, \
                                lon0, lat0, z0, zpar, do3d, doturb, name, N=N, \
                                grid=grid, dostream=dostream, T0=T0, U=U, V=V)

        # # nan out drifters once they have been used to run tracpy
        # startlon[idrifts] = np.nan
        # startlat[idrifts] = np.nan
        # startdate[idrifts] = np.nan

        # pdb.set_trace()

if __name__ == "__main__":
    run()    