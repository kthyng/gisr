'''
Functions to initialize various numerical experiments.

Make a new init_* for your application.

loc     Path to directory of grid and output files
nsteps  Number of steps to do between model outputs (iter in tracmass)
ndays   number of days to track the particles from start date
ff      ff=1 to go forward in time and ff=-1 for backward in time
date    Start date in datetime object
tseas   Time between outputs in seconds
ah      Horizontal diffusion in m^2/s. 
        See project values of 350, 100, 0, 2000. For -turb,-diffusion
av      Vertical diffusion in m^2/s.
do3d    for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
doturb  turbulence/diffusion flag. 
        doturb=0 means no turb/diffusion,
        doturb=1 means adding parameterized turbulence
        doturb=2 means adding diffusion on a circle
        doturb=3 means adding diffusion on an ellipse (anisodiffusion)
lon0    Drifter starting locations in x/zonal direction.
lat0    Drifter starting locations in y/meridional direction.
z0/zpar Then z0 should be an array of initial drifter depths. 
        The array should be the same size as lon0 and be negative
        for under water. Currently drifter depths need to be above 
        the seabed for every x,y particle location for the script to run.
        To do 3D but start at surface, use z0=zeros(ia.shape) and have
         either zpar='fromMSL'
        choose fromMSL to have z0 starting depths be for that depth below the base 
        time-independent sea level (or mean sea level).
        choose 'fromZeta' to have z0 starting depths be for that depth below the
        time-dependent sea surface. Haven't quite finished the 'fromZeta' case.
        Then: 
        set z0 to 's' for 2D along a terrain-following slice
         and zpar to be the index of s level you want to use (0 to km-1)
        set z0 to 'rho' for 2D along a density surface
         and zpar to be the density value you want to use
         Can do the same thing with salinity ('salt') or temperature ('temp')
         The model output doesn't currently have density though.
        set z0 to 'z' for 2D along a depth slice
         and zpar to be the constant (negative) depth value you want to use
        To simulate drifters at the surface, set z0 to 's' 
         and zpar = grid['km']-1 to put them in the upper s level
         z0='s' is currently not working correctly!!!
         In the meantime, do surface using the 3d set up option but with 2d flag set
xp      x-locations in x,y coordinates for drifters
yp      y-locations in x,y coordinates for drifters
zp      z-locations (depths from mean sea level) for drifters
t       time for drifter tracks
name    Name of simulation to be used for netcdf file containing final tracks

'''

import numpy as np
import os
import netCDF4 as netCDF
import pdb
import glob
from datetime import datetime, timedelta
from matplotlib.mlab import *
import inout
import tools

units = 'seconds since 1970-01-01'

def sensitivity(loc=None, nsteps=None, ff=None, ah=None, grid=None, nlon=None, nlat=None, doturb=None, name=None):
    '''
    A drifter test using TXLA model output. 
    The comparison case for this simulation is 2D (do3d=0) 
    with no turbulence/diffusion (doturb=0).
    Drifters are started at the surface and run forward
    for ten days (ndays=10) from 11/25/09 (in date). Compare results with figure in examples/test1.png.

    Optional inputs for making tests easy to run:
        loc             'thredds' or 'local', default = 'thredds'
        nsteps          Number of particle steps to record between model outputs
                        Default = 5
        ff              Backward (-1) or forward (1) in time. Default is forward (1).
        ah              Horizontal viscosity, default = 5
        grid            If input, will not redo this step. Default is to load in grid.
        nlon, nlat      Number of drifters to use in the lon/lat direction in seed array
                        Default = 110, 98 (10 km spacing)
        doturb          What, if any, subgrid parameterization to use. Default is 'none'
        name            Specific name for track and figure files. Default is 'temp'
    '''

    # Location of TXLA model output
    # file and then grid. 
    # 0150 file goes from (2009, 11, 19, 12, 0) to (2009, 12, 6, 0, 0)
    if loc is None or loc == 'thredds':
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'      
    elif loc is 'local':
    # Location of TXLA model output
        if 'rainier' in os.uname():
            loc = '/Users/kthyng/Documents/research/postdoc/' # for model outputs
        elif 'hafen.tamu.edu' in os.uname():
            loc = '/home/kthyng/shelf/' # for model outputs

    # Initialize parameters
    if nsteps is None:
        nsteps = 5
    else:
        nsteps = nsteps

    ndays = 16
    if ff is None:
        ff = 1
    else:
        ff = ff
    # Start date
    # date = datetime(2009, 11, 19, 12)
    date = datetime(2009, 11, 20, 0)

    # Time between outputs
    # Dt = 14400. # in seconds (4 hours), nc.variables['dt'][:] 
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    if ah is None:
        ah = 5. #100.
    else:
        ah = ah
    
    av = 1.e-5 # m^2/s, or try 5e-6

    # grid = netCDF.Dataset(loc+'grid.nc')
    # lonr = grid.variables['lon_rho'][:]
    # latr = grid.variables['lat_rho'][:]
    if grid is None:
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    if nlon is None:
        nlon = 110
    else:
        nlon = nlon
    if nlat is None:
        nlat = 98
    else:
        nlat = nlat
    lon0,lat0 = np.meshgrid(np.linspace(-98.5,-87.5,nlon),np.linspace(22.5,31,nlat)) # whole domain, 10 km

    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tools.check_points(lon0,lat0,grid)

    ## Choose method for vertical placement of drifters
    z0 = 's'  #'salt' #'s' #'z' #'salt' #'s' 
    zpar = 29 #30 #29 #-10 #grid['km']-1 # 30 #grid['km']-1
    # Do the following two for a 3d simulation
    # z0 = np.ones(xstart0.shape)*-40 #  below the surface
    # zpar = 'fromMSL' 
    # pdb.set_trace()

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    # turbulence/diffusion flag. doturb=0 means no turb/diffusion,
    # doturb=1 means adding parameterized turbulence
    # doturb=2 means adding diffusion on a circle
    # doturb=3 means adding diffusion on an ellipse (anisodiffusion)
    if doturb is None:
        doturb = 0
    else:
        doturb = doturb

    # simulation name, used for saving results into netcdf file
    if name is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = name

    return loc,nsteps,ndays,ff,date,tseas,ah,av,lon0,lat0,z0,zpar,do3d,doturb,name,grid

def outer_f(date=None, grid=None):
    '''
    Initialization for seeding drifters along the outside of the 
    numerical domain to be run forward.

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0). If date not input,
                name will be 'temp' 
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ndays = 60
    ff = 1 # This is a backward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    # Assign them in grid coordinates
    i0a = np.arange(0,grid['imt'],6)
    j0a = np.ones(i0a.shape)*3
    j0b = np.arange(0,grid['jmt'],6)
    i0b = np.ones(j0b.shape)*3
    # i0c = np.arange(0,grid['imt'],6)
    # j0c = np.ones(i0a.shape)*(grid['jmt']-3)
    j0d = np.arange(0,grid['jmt'],6)
    i0d = np.ones(j0d.shape)*(grid['imt']-3)
    i0 = np.hstack((i0a,i0b,i0d))
    j0 = np.hstack((j0a,j0b,j0d))

    # Convert to lon/lat
    lon0, lat0, _ = tools.interpolate2d(i0, j0, grid,'m_ij2ll',
                                        mode='constant', cval=np.nan)

    # pdb.set_trace()
    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tools.check_points(lon0,lat0,grid)
    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'outer_f/' + date.isoformat()[0:10] 

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid

def galv_b(date=None, grid=None):
    '''
    Initialization for seeding drifters near Galveston Bay to be run
    backward.

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0). If date not input,
                name will be 'temp' 
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ndays = 90
    ff = -1 # This is a backward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    lon0,lat0 = np.meshgrid(np.linspace(-95.3,-94.3,30), 
                            np.linspace(28.6,29.6,30))

    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tools.check_points(lon0,lat0,grid)

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 1
    # convert date to number
    datenum = netCDF.date2num(date, units)
    # Number of model outputs to use
    tout = np.int((ndays*(24*3600))/tseas)
    # Figure out what files will be used for this tracking - to get tinds for
    # the following calculation
    nc, tinds = inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()
    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = (abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/N
    # Initialize arrays of lon0, lat0 and U, V for full number of drifters
    lon0 = np.ones(N,order='F')*lon0
    lat0 = np.ones(N,order='F')*lat0
    T0 = np.ones(N,order='F')*T0

    # Initialize the arrays to save the transports on the grid in the loop.
    # These arrays aggregate volume transport when a drifter enters or exits a grid cell
    # These should start at zero since we don't know which way things will travel yet
    U = np.ma.zeros(grid['xu'].shape,order='F')
    V = np.ma.zeros(grid['xv'].shape,order='F')

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'galv_b/' + date.isoformat()[0:13] 

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid, dostream, T0, U, V

def galv_fromb2f(date=None, grid=None):
    '''
    Initialization for seeding drifters from backward final locations 
    near Galveston Bay to be run
    forward.

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0). If date not input,
                name will be 'temp' 
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ndays = 90
    ff = 1 # This is a backward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'galv_fromb2f/' + date.isoformat()[0:13] 

    return loc, nsteps, ndays, ff, tseas, ah, av, \
            z0, zpar, do3d, doturb, name, grid

def galv_f(date=None, grid=None):
    '''
    Initialization for seeding drifters near Galveston Bay to be run
    forward.

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0). If date not input,
                name will be 'temp' 
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ndays = 60 # run drifters for ndays
    ff = 1 # This is a backward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    lon0,lat0 = np.meshgrid(np.linspace(-95.3,-94.3,30), 
                            np.linspace(28.6,29.6,30))

    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tools.check_points(lon0,lat0,grid)

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 1
    # convert date to number
    datenum = netCDF.date2num(date, units)
    # Number of model outputs to use
    tout = np.int((ndays*(24*3600))/tseas)
    # Figure out what files will be used for this tracking - to get tinds for
    # the following calculation
    nc, tinds = inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()
    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = (abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/N
    # Initialize arrays of lon0, lat0 and U, V for full number of drifters
    lon0 = np.ones(N,order='F')*lon0
    lat0 = np.ones(N,order='F')*lat0
    T0 = np.ones(N,order='F')*T0

    # Initialize the arrays to save the transports on the grid in the loop.
    # These arrays aggregate volume transport when a drifter enters or exits a grid cell
    # These should start at zero since we don't know which way things will travel yet
    U = np.ma.zeros(grid['xu'].shape,order='F')
    V = np.ma.zeros(grid['xv'].shape,order='F')

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'galv_f/' + date.isoformat()[0:13] 

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid, dostream, T0, U, V

def bara_b(ndatum=0, hour=0, grid=None):
    '''
    Initialization for seeding drifters near Barataria Bay to be run
    backward. Start a 2d gaussian of drifters in x y backward started
    over a range of 2 days with the number of drifters run selected 
    using a time gaussian.

    Optional inputs for making tests easy to run:
        ndatum  Which oil datum is being run, in order from 0 to 55.
                Default is 0.
        hour    What hour of 48 this simulation is being run for.
                Default is 0.
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    import bara_data

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Get lon/lat/time info from bara_data.py
    lon, lat, date = bara_data.retrieve()
    # just use the lon/lat/date for this oil datum
    lon = lon[ndatum]
    lat = lat[ndatum]
    # run for a variety of times before oil was found
    date = netCDF.num2date(date[ndatum],units) - timedelta(hours=hour)

    # when oil started spilling, 4/20/10 9:45PM CDT to UTC
    enddate = datetime(2010, 4, 20, 9+12, 45, 0) + timedelta(hours=5)

    # Number of days to run depends on date of data relative to spill start
    ndays = (date - enddate).total_seconds()/(3600.*24)

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ff = -1 # This is a backward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 100.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # Initialize drifters in space using 2d gaussian.
    # Select out the lon/lat for test
    dlon = 0.1; dlat = 0.14 # delta degree distances for starting particles
    # Time Gaussian to set number of drifters used in (x,y)
    H = np.arange(48) # hours in 2 days
    mu = 24 # 1 day into the 2 days of simulation starts is the mean
    sigma = 16. # Standard deviation
    # pdb.set_trace()
    N = 30*np.exp(-(H-mu)**2/(2*sigma**2))
    # Choose N value for hour
    Nh = np.floor(N[H==hour])

    lon0, lat0 = tools.seed(lon, lat, dlon=dlon, dlat=dlat, N=Nh)

    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tools.check_points(lon0,lat0,grid)

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 1

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'bara_b/turb/' + date.isoformat()

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid

def dwh_f(date=None, grid=None):
    '''
    Initialization for seeding drifters near the Deepwater Horizon
    accident site to be run forward.

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0). If date not input,
                name will be 'temp' 
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ndays = 90
    ff = 1 # This is a forward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    lon0,lat0 = np.meshgrid(np.linspace(-88.81-.4,-88.81+.4,15), 
                            np.linspace(28.1,28.9,15))

    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tools.check_points(lon0,lat0,grid)

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'dwh_f/' + date.isoformat()[0:10] 

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid

def dwh_stream_f(date, N, grid=None):
    '''
    Initialization for seeding drifters near the Deepwater Horizon
    accident site to be run forward. for lagrangian streamfunctions

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0).
        N       Number of drifters to seed
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ndays = 90
    ff = 1 # This is a forward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 20.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    # These are the actual blowout coords for DWH (shifted to be in domain)
    lon0 = np.array([-88.36594444444444-.15])
    lat0 = np.array([28.73813888888889+.15])
    # # Want to use the coords for the nearby east/north grid cell walls for
    # # the proper transport calculation. Index: (617, 6)
    # lon0 = np.array([-88.50377451891225])
    # lat0 = np.array([28.888202179642793])
    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]
    # Change to get positions at the center of the given cell
    lon0, lat0, _ = tools.interpolate2d(ia - 0.5, ja - 0.5, grid, 'm_ij2ll')

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 1

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 1
    # convert date to number
    datenum = netCDF.date2num(date, units)
    # Number of model outputs to use
    tout = np.int((ndays*(24*3600))/tseas)
    # Figure out what files will be used for this tracking - to get tinds for
    # the following calculation
    nc, tinds = inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()
    # Save initial volume transport of each drifter. Initial volume is equal to
    # the velocity of the initial drifter locations times the flux and divided
    # by the number of drifters in that cell.
    # All drifters in this case are initialized in the same grid cell. T is overall
    # volume transport.
    # pdb.set_trace()
    # if uf[ia, ja, 0]>0: # if positive u flux, start u drifters on positive side of cell
    #     iau = ia
    #     jau = ja
    # else: # if negative u flux, start u drifters on negative side of cell
    #     iau = ia-1
    #     jau = ja
    # lon0u, lat0u = grid['basemap'](grid['xpsi'][iau, jau], grid['ypsi'][iau, jau], \
    #                                 inverse=True)
    # if vf[ia, ja, 0]>0: # if positive v flux, start v drifters on positive side of cell
    #     iav = ia
    #     jav = ja
    # else: # if negative v flux, start v drifters on negative side of cell
    #     iav = ia
    #     jav = ja-1
    # lon0v, lat0v = grid['basemap'](grid['xpsi'][iav, jav], grid['ypsi'][iav, jav], \
    #                                 inverse=True)
    # # Base the number of drifters used in each direction on the ratio of 
    # # transport in each direction
    # Nu = int(N*(abs(uf[iau, jau, 0])/(abs(uf[iau, jau, 0])+abs(vf[iav, jav, 0]))))
    # Nv = int(N*(abs(vf[iav, jav, 0])/(abs(uf[iau, jau, 0])+abs(vf[iav, jav, 0]))))
    # # Initial volume transport for u and v-moving drifters
    # U0 = uf[ia, ja, 0]/N
    # V0 = vf[ia, ja, 0]/N
    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = (abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/N
    # Initialize arrays of lon0, lat0 and U, V for full number of drifters
    lon0 = np.ones(N,order='F')*lon0
    lat0 = np.ones(N,order='F')*lat0
    T0 = np.ones(N,order='F')*T0
    # U0 = np.ones(N,order='F')*U0
    # V0 = np.ones(N,order='F')*V0
    # # Initialize arrays of lon0, lat0 and U, V for full number of drifters
    # lon0 = np.concatenate((np.ones(Nu)*lon0u, np.ones(Nv)*lon0v))
    # lat0 = np.concatenate((np.ones(Nu)*lat0u, np.ones(Nv)*lat0v))
    # U0 = np.concatenate((np.ones(Nu)*U0, np.ones(Nv)*V0))

    # Initialize the arrays to save the transports on the grid in the loop.
    # These arrays aggregate volume transport when a drifter enters or exits a grid cell
    # These should start at zero since we don't know which way things will travel yet
    U = np.ma.zeros(grid['xu'].shape,order='F')
    # Urho[ia, ja] = Urho[ia, ja] + U0.sum()
    V = np.ma.zeros(grid['xv'].shape,order='F')
    # Vrho[ia, ja] = Vrho[ia, ja] + V0.sum()

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'dwh_stream_f/' + date.isoformat()[0:13] + 'N' + str(N)

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid, dostream, T0, U, V

def bara_stream_b(date, runend, N, grid=None):
    '''
    Initialization for seeding drifters near the Deepwater Horizon
    accident site to be run forward. for lagrangian streamfunctions

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0).
        N       Number of drifters to seed
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ff = -1 # This is a forward-moving simulation

    # Number of days to run depends on date of data relative to spill start
    ndays = (date - runend).total_seconds()/(3600.*24)

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 20.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    # Outside Barataria Bay
    lon0 = np.array([-89.88])
    lat0 = np.array([29.22])
    # # Want to use the coords for the nearby east/north grid cell walls for
    # # the proper transport calculation. Index: (617, 6)
    # lon0 = np.array([-88.50377451891225])
    # lat0 = np.array([28.888202179642793])
    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]
    # Change to get positions at the center of the given cell
    lon0, lat0, _ = tools.interpolate2d(ia - 0.5, ja - 0.5, grid, 'm_ij2ll')

    # surface drifters
    z0 = 's'  
    zpar = 29 

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 1

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 1
    # convert date to number
    datenum = netCDF.date2num(date, units)
    # Number of model outputs to use
    tout = np.int((ndays*(24*3600))/tseas)
    # Figure out what files will be used for this tracking - to get tinds for
    # the following calculation
    nc, tinds = inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = inout.readfields(tinds[0],grid,nc,z0,zpar)
    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = (abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/N
    # Initialize arrays of lon0, lat0 and U, V for full number of drifters
    lon0 = np.ones(N,order='F')*lon0
    lat0 = np.ones(N,order='F')*lat0
    T0 = np.ones(N,order='F')*T0

    # Initialize the arrays to save the transports on the grid in the loop.
    # These arrays aggregate volume transport when a drifter enters or exits a grid cell
    U = np.zeros(grid['xu'].shape,order='F')
    V = np.zeros(grid['xv'].shape,order='F')

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'bara_stream_b/' + date.isoformat()[0:13] + 'N' + str(N)

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid, dostream, T0, U, V
