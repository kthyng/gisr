'''
The functions here calculate and save netcdf files of the velocity field 
for a gyre test of tracpy. 
'''

from __future__ import division
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

def make_fields(grid):
    '''
    Make u,v fields for gyre test
    Input:
        grid    as read in from tracpy.inout.readgrid()
    '''

    # Center of gyre in lon/lat
    lon0 = -93.5; lat0 = 28.5

    # Convert to x/y
    x0, y0 = grid['basemap'](lon0, lat0)

    # Size of gyre
    L = 30000 # meters

    # Velocity fields, on staggered grid
    xu = grid['xu']
    u = (grid['yu']-y0)/L**2*np.exp(-1./(2*L**2)*((grid['xu']-x0)**2 + (grid['yu']-y0)**2))
    v = -(grid['xv']-x0)/L**2*np.exp(-1./(2*L**2)*((grid['xv']-x0)**2 + (grid['yv']-y0)**2))

    # Also need time
    t = np.ones(24*5/4)*np.nan
    for i,hour in enumerate(np.arange(0, 24*5, 4)): # five days, 4 hours at a time
        t[i] = netCDF.date2num(datetime(2012, 5, 23, 0, 0) + timedelta(hours=int(hour)),units)
    # pdb.set_trace()
    # Make u and v have a time dimension
    u = u.reshape((1,u.shape[0],u.shape[1])).repeat(t.size,axis=0)
    v = v.reshape((1,v.shape[0],v.shape[1])).repeat(t.size,axis=0)

    return u, v, t


def write_files(uin, vin, tin):

    # pdb.set_trace()

    # Open file for writing.
    rootgrp = netCDF.Dataset('projects/gyre/ocean_his_0001.nc', 'w', format='NETCDF3_64BIT')

    # find dimensions
    el = vin.shape[1]
    xl = uin.shape[2]
    nt = tin.size

    # Define dimensions
    rootgrp.createDimension('el', el) # length of array in y direction
    rootgrp.createDimension('xl', xl) # x direction
    rootgrp.createDimension('elm1', el-1) # length of array in y direction
    rootgrp.createDimension('xlm1', xl-1) # x direction
    rootgrp.createDimension('nt', nt) # time length

    # Create variable
    u = rootgrp.createVariable('u','f8',('nt','elm1','xl')) # 64-bit floating point
    v = rootgrp.createVariable('v','f8',('nt','el','xlm1')) # 64-bit floating point
    ocean_time = rootgrp.createVariable('ocean_time','f8',('nt'))

    # # Set some attributes
    # u.long_name = 'longitudinal position of drifter'
    # lonp.units = 'degrees'
    # lonp.time = 'tp'

    # Write data to netCDF variables
    u[:] = uin
    v[:] = vin
    ocean_time[:] = tin
    rootgrp.close()

def uv(grid):

    # If the netCDF doesn't exist, make the file
    # if not os.path.exists('projects/gyre'):
    #     os.makedirs('projects/gyre')
    # if not os.path.exists('projects/gyre/ocean_his_0001.nc'):
    u, v, t = make_fields(grid)
    write_files(u, v, t)

def init():

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ff = 1
    ndays = 5

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    lon0 = np.linspace(-93.5,-93,10)
    lat0 = np.ones(lon0.shape)*28.5

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0
    dostream = 0

    name = 'gyre'

    return nsteps, ndays, ff, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, dostream, name


def run():

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Read in grid
    grid = tracpy.inout.readgrid(loc)

    # Create gyre field velocities if don't exist yet
    uv(grid)

    # Start from the beginning and add days on for loop
    date = datetime(2012, 5, 23, 0, 0)

    # Read in simulation initialization
    nsteps, ndays, ff, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, dostream, name = init()

    # Run tracpy
    loc = 'projects/gyre/' # loc for model output
    lonp, latp, zp, t, grid \
        = tracpy.run.run(loc, nsteps, ndays, ff, date, tseas, ah, av, \
                            lon0, lat0, z0, zpar, do3d, doturb, name, \
                            grid=grid, dostream=dostream)



if __name__ == "__main__":
    run()    