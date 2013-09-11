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
import octant
from matplotlib.mlab import find
from mpl_toolkits.basemap import Basemap

units = 'seconds since 1970-01-01'

def make_grid():
    
    # Make grid (did this on hafen)
    lon, lat = mgrid[-95.:-94.0:.001, 29:30:.001]
    llcrnrlon=-95; llcrnrlat=29; 
    urcrnrlon=-94; urcrnrlat=30; projection='lcc'
    lat_0=29.5; lon_0=-94.5; resolution='i'; area_thresh=0.
    basemap = Basemap(llcrnrlon=llcrnrlon,
                 llcrnrlat=llcrnrlat,
                 urcrnrlon=urcrnrlon,
                 urcrnrlat=urcrnrlat,
                 projection=projection,
                 lat_0=lat_0,
                 lon_0=lon_0,
                 resolution=resolution,
                 area_thresh=area_thresh)
    grd = octant.grid.CGrid_geo(lon, lat, basemap)
    grd.h=np.ones(grd.lat_rho.shape)*100
    # cd to octant directory
    # cd /usr/local/lib/python2.7/lib/python2.7/site-packages/octant
    # import roms
    roms.write_grd(grd,filename='/home/kthyng/projects/gyre/grid.nc',verbose=True)

    s_win = np.array([-1.,0.])
    Cs_win = np.array([-1.,0.])
    hcin = 0.
    theta_sin = 1e-4
    theta_bin = 1.

    # Append to grid
    # Open file for writing.
    rootgrp = netCDF.Dataset('/home/kthyng/projects/gyre/grid.nc', 'a', format='NETCDF3_64BIT')

    # find dimensions
    zl = 1

    # Define dimensions
    rootgrp.createDimension('zlp1', zl+1) # z direction

    # Create variable
    s_w = rootgrp.createVariable('s_w','f8',('zlp1'))
    Cs_w = rootgrp.createVariable('Cs_w','f8',('zlp1'))
    hc = rootgrp.createVariable('hc','f8')
    theta_s = rootgrp.createVariable('theta_s','f8')
    theta_b = rootgrp.createVariable('theta_b','f8')

    # Write data to netCDF variables
    s_w[:] = s_win
    Cs_w[:] = Cs_win
    hc[:] = hcin
    theta_s[:] = theta_sin
    theta_b[:] = theta_bin
    rootgrp.close()


# llcrnrlon=-95; llcrnrlat=29; 
# urcrnrlon=-94; urcrnrlat=30; projection='lcc'
# lat_0=29.5; lon_0=-94.5; resolution='i'; area_thresh=0.
# basemap = Basemap(llcrnrlon=llcrnrlon,
#              llcrnrlat=llcrnrlat,
#              urcrnrlon=urcrnrlon,
#              urcrnrlat=urcrnrlat,
#              projection=projection,
#              lat_0=lat_0,
#              lon_0=lon_0,
#              resolution=resolution,
#              area_thresh=area_thresh)

# # grdll = octant.grid.CGrid_geo(grd,x,y, basemap)
# lon = (llcrnrlon,llcrnrlon,urcrnrlon,urcrnrlon)
# lat = (llcrnrlat,llcrnrlat,urcrnrlat,urcrnrlat)
# beta = [1.,1.,1.,1.]
# grd = octant.grid.Gridgen(lon, lat, beta, (32,32), proj=basemap)


# proj = Basemap(projection='lcc',
#                resolution='i',
#                llcrnrlon=-72.0,
#                llcrnrlat= 40.0,
#                urcrnrlon=-63.0,
#                urcrnrlat=47.0,
#                lat_0=43.0,
#                lon_0=-62.5)

# lon = (-71.977385177601761, -70.19173825913137,
#        -63.045075098584945,-64.70104074097425)
# lat = (42.88215610827428, 41.056141745853786,
#        44.456701607935841, 46.271758064353897)
# beta = [1.0, 1.0, 1.0, 1.0]

# grd = octant.grid.Gridgen(lon, lat, beta, (32, 32), proj=proj)

def make_fields(grid, ndays):
    '''
    Make u,v fields for gyre test
    Input:
        grid    as read in from tracpy.inout.readgrid()
    '''
    # pdb.set_trace()
    # Center of gyre in lon/lat
    lon0 = -94.5; lat0 = 29.5

    # Convert to x/y
    x0, y0 = grid['basemap'](lon0, lat0)

    # Radius of gyre
    L = 10000 # meters

    # Strength of gyre
    psi = 1000000

    # Velocity fields, on staggered grid
    u = psi*(grid['yu'].T-y0)/L**2*np.exp(-1./(2*L**2)*((grid['xu'].T-x0)**2 + (grid['yu'].T-y0)**2))
    v = -psi*(grid['xv'].T-x0)/L**2*np.exp(-1./(2*L**2)*((grid['xv'].T-x0)**2 + (grid['yv'].T-y0)**2))

    # Also need time
    t = np.ones(24*ndays/4+1)*np.nan
    for i,hour in enumerate(np.arange(0, 24*ndays+1, 4)): # five days, 4 hours at a time
        t[i] = netCDF.date2num(datetime(2012, 5, 23, 0, 0) + timedelta(hours=int(hour)),units)
    # pdb.set_trace()
    # Make u and v have a time dimension
    u = u.reshape((1,u.shape[0],u.shape[1])).repeat(t.size,axis=0)
    v = v.reshape((1,v.shape[0],v.shape[1])).repeat(t.size,axis=0)

    # also need depth
    # Make u and v have a k dimension
    u = u.reshape((u.shape[0],1,u.shape[1],u.shape[2]))
    v = v.reshape((u.shape[0],1,v.shape[1],v.shape[2]))

    # also need zeta
    zeta = np.zeros((u.shape[0],u.shape[2],v.shape[3]))

    # other grid parameters:
    s_w = np.array([0,1])
    Cs_w = np.array([0,1])
    hc = 0.
    theta_s = 0.
    theta_b = 1e-4

    return u, v, t, zeta, s_w, Cs_w, hc, theta_s, theta_b


def write_files(uin, vin, tin, zetain, s_win, Cs_win, hcin, 
                theta_sin, theta_bin):

    # pdb.set_trace()

    # Open file for writing.
    rootgrp = netCDF.Dataset('projects/gyre/ocean_his_0001.nc', 'w', format='NETCDF3_64BIT')

    # find dimensions
    el = uin.shape[2]
    xl = vin.shape[3]
    zl = uin.shape[1]
    nt = tin.size

    # Define dimensions
    rootgrp.createDimension('el', el) # length of array in y direction
    rootgrp.createDimension('xl', xl) # x direction
    rootgrp.createDimension('zl', zl) # z direction
    rootgrp.createDimension('zlp1', zl+1) # z direction
    rootgrp.createDimension('elm1', el-1) # length of array in y direction
    rootgrp.createDimension('xlm1', xl-1) # x direction
    rootgrp.createDimension('nt', None) # time length, setting this as None makes
    # this unlimited and also makes it the default aggdim for MFDataset

    # Create variable
    u = rootgrp.createVariable('u','f8',('nt','zl','el','xlm1')) # 64-bit floating point
    v = rootgrp.createVariable('v','f8',('nt','zl','elm1','xl')) # 64-bit floating point
    zeta = rootgrp.createVariable('zeta','f8',('nt','el','xl')) # 64-bit floating point
    ocean_time = rootgrp.createVariable('ocean_time','f8',('nt'))
    s_w = rootgrp.createVariable('s_w','f8',('zlp1'))
    Cs_w = rootgrp.createVariable('Cs_w','f8',('zlp1'))
    hc = rootgrp.createVariable('hc','f8')
    theta_s = rootgrp.createVariable('theta_s','f8')
    theta_b = rootgrp.createVariable('theta_b','f8')

    # # Set some attributes
    # u.long_name = 'longitudinal position of drifter'
    # lonp.units = 'degrees'
    # lonp.time = 'tp'

    # Write data to netCDF variables
    u[:] = uin
    v[:] = vin
    zeta[:] = zetain
    ocean_time[:] = tin
    s_w[:] = s_win
    Cs_w[:] = Cs_win
    hc[:] = hcin
    theta_s[:] = theta_sin
    theta_b[:] = theta_bin
    rootgrp.close()

def uv(grid, ndays):

    # If the netCDF doesn't exist, make the file
    # if not os.path.exists('projects/gyre'):
    #     os.makedirs('projects/gyre')
    # if not os.path.exists('projects/gyre/ocean_his_0001.nc'):
    u, v, t, zeta, s_w, Cs_w, hc, theta_s, theta_b = make_fields(grid, ndays)
    write_files(u, v, t, zeta, s_w, Cs_w, hc, theta_s, theta_b)

def init():

    # Initialize parameters
    nsteps = 20 # time interpolation steps
    ff = 1
    ndays = 10

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    lon0 = np.linspace(-94.5,-94,10)
    lat0 = np.ones(lon0.shape)*29.5

    # surface drifters
    z0 = 's'  
    zpar = 0

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0
    dostream = 0

    name = 'gyre/gyre'

    return nsteps, ndays, ff, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, dostream, name


def run():

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/gyre'):
        os.makedirs('tracks/gyre')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/gyre'):
        os.makedirs('figures/gyre')

    # Location of TXLA model output
    loc = 'projects/gyre/'# 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Read in grid
    grid = tracpy.inout.readgrid(loc)

    # Read in simulation initialization
    nsteps, ndays, ff, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, dostream, name = init()

    # Create gyre field velocities if don't exist yet
    uv(grid, ndays)

    # Start from the beginning and add days on for loop
    date = datetime(2012, 5, 23, 0, 0)
    # pdb.set_trace()
    # Run tracpy
    loc = 'projects/gyre/' # loc for model output
    lonp, latp, zp, t, grid \
        = tracpy.run.run(loc, nsteps, ndays, ff, date, tseas, ah, av, \
                            lon0, lat0, z0, zpar, do3d, doturb, name, \
                            grid=grid, dostream=dostream)

    # Read in and plot tracks
    d = netCDF.Dataset('tracks/' + name + '.nc')
    lonp = d.variables['lonp'][:]
    latp = d.variables['latp'][:]
    tracpy.plotting.tracks(lonp, latp, name, grid=grid)
    tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
    d.close()


if __name__ == "__main__":
    run()    