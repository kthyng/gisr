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
import tracpy

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
        grid = tracpy.inout.readgrid(loc)
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
    lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)

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
        grid = tracpy.inout.readgrid(loc)
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
    lon0, lat0, _ = tracpy.tools.interpolate2d(i0, j0, grid,'m_ij2ll',
                                        mode='constant', cval=np.nan)

    # pdb.set_trace()
    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)
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
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    # lon0,lat0 = np.meshgrid(np.linspace(-95.3,-94.3,30), 
    #                         np.linspace(28.6,29.6,30))
    # These are from the first simulation of Matt's, directly (5-23-10-T00)
    lon0 = np.array([-95.0241399 , -94.98965578, -94.95517325, -94.92069183,
       -94.88620883, -94.85172509, -94.81724335, -94.78275982,
       -94.7482779 , -94.71379507, -94.67931158, -95.16206948,
       -95.12758717, -95.09310538, -95.05862145, -95.02413794,
       -94.98965701, -94.95513931, -94.9206318 , -94.88614683,
       -94.85167001, -94.81717435, -94.78272419, -94.7482423 ,
       -94.71379539, -94.67931198, -94.64482931, -94.61034688,
       -95.2655187 , -95.2310361 , -95.19655087, -95.16207005,
       -95.12758782, -95.09305431, -95.05857084, -95.02410113,
       -94.98962946, -94.95511785, -94.92061517, -94.88611981,
       -94.85162388, -94.81713473, -94.78264482, -94.7481572 ,
       -94.71366432, -94.67922921, -94.64482815, -94.61034697,
       -94.57586365, -95.3000017 , -95.2655181 , -95.2310333 ,
       -95.19653484, -95.1620388 , -95.12755681, -95.0930927 ,
       -95.0586025 , -95.02410465, -94.98961709, -94.95512324,
       -94.92061694, -94.88611758, -94.85162163, -94.81713356,
       -94.78263082, -94.74814088, -94.71364736, -94.67914544,
       -94.64466756, -94.61020375, -94.57586434, -94.5413805 ,
       -94.50689871, -95.29999722, -95.26551799, -95.23105674,
       -95.19657054, -95.16205672, -95.12757422, -95.09306828,
       -95.05860208, -95.02411644, -94.9896116 , -94.95511559,
       -94.92060986, -94.8861118 , -94.85162132, -94.81713422,
       -94.78264541, -94.74815372, -94.71370889, -94.67918341,
       -94.64467528, -94.61020068, -94.57572683, -94.54138145,
       -94.50689755, -94.47241587, -94.4379322 , -95.26555061,
       -95.23103205, -95.19653341, -95.16204189, -95.12755044,
       -95.0930626 , -95.0585794 , -95.02409261, -94.98960462,
       -94.95510265, -94.9206037 , -94.88610817, -94.85162115,
       -94.81713221, -94.78264374, -94.74815832, -94.71366663,
       -94.67925681, -94.64469702, -94.61022039, -94.57573748,
       -94.54125526, -94.50679148, -94.47241496, -94.43793323,
       -94.40344925, -94.36896757, -95.19649488, -95.1620103 ,
       -95.12752738, -95.0930437 , -95.05855748, -95.02407384,
       -94.98957927, -94.95508959, -94.92060034, -94.88611059,
       -94.85161738, -94.81712686, -94.78264129, -94.74815407,
       -94.71366744, -94.67918942, -94.64470334, -94.61022465,
       -94.57583904, -94.54130193, -94.50677423, -94.47231023,
       -94.43793233, -94.40345044, -94.36896651, -95.16199311,
       -95.12750233, -95.09302201, -95.05853444, -95.02404783,
       -94.98955697, -94.95507538, -94.9205921 , -94.88610191,
       -94.85161473, -94.81712553, -94.78263937, -94.74815334,
       -94.71366913, -94.67929834, -94.64471377, -94.61023533,
       -94.57575945, -94.54126356, -94.5067751 , -94.47228261,
       -94.43793208, -94.40344962, -94.36896753, -95.1275201 ,
       -95.09300681, -95.05853101, -95.0240449 , -94.989561  ,
       -94.95507202, -94.92058516, -94.88609511, -94.85161244,
       -94.8171244 , -94.78263975, -94.74815675, -94.7136751 ,
       -94.67919383, -94.64470942, -94.61022975, -94.57574723,
       -94.54126392, -94.50677672, -94.47228188, -94.43779522,
       -94.40344873, -94.36896673, -95.09307524, -95.0585704 ,
       -95.02404034, -94.98955968, -94.95506705, -94.9205787 ,
       -94.88609684, -94.85161301, -94.81712623, -94.78264071,
       -94.74816118, -94.71367794, -94.67919103, -94.64470725,
       -94.61022116, -94.57574044, -94.54125494, -94.50676922,
       -94.4722896 , -94.43779573, -94.40330904, -94.36896617,
       -95.0240494 , -94.98955984, -94.95506912, -94.92058844,
       -94.88610242, -94.8516158 , -94.81713094, -94.78264925,
       -94.74816199, -94.71367295, -94.67919152, -94.64470357,
       -94.61021991, -94.57573713, -94.54125231, -94.50676919,
       -94.47227957, -94.43779669, -94.40330223, -94.36896737,
       -94.33448354, -94.98957868, -94.955084  , -94.92059823,
       -94.88611136, -94.85162801, -94.81714058, -94.78265147,
       -94.74816622, -94.71367877, -94.67919459, -94.64470943,
       -94.6102225 , -94.57573871, -94.54125596, -94.50676809,
       -94.47227971, -94.43779329, -94.40330804, -94.36896646,
       -94.33448473, -94.920616  , -94.88612119, -94.85163229,
       -94.81714512, -94.78265761, -94.74816964, -94.71379487,
       -94.67919069, -94.64470996, -94.61022605, -94.5757423 ,
       -94.54125522, -94.50677662, -94.47228688, -94.43779977,
       -94.40331201, -94.36896743, -94.33448378, -94.88620814,
       -94.85165828, -94.81716063, -94.78267083, -94.74824137,
       -94.71368664, -94.67920442, -94.64472177, -94.61022466,
       -94.57574458, -94.54126575, -94.50678072, -94.47229048,
       -94.43780518, -94.40331805, -94.36896525, -94.3344843 ,
       -94.8171858 , -94.78268892, -94.74820967, -94.71371169,
       -94.67922681, -94.64471383, -94.61023142, -94.57575611,
       -94.54126891, -94.50678384, -94.47229335, -94.43781414,
       -94.40331992, -94.36882327, -94.33448285, -94.74824019,
       -94.71373235, -94.67923883, -94.64476681, -94.61023978,
       -94.57575757, -94.54128176, -94.50679793, -94.47231165,
       -94.43781791, -94.40331198, -94.36896724, -94.33448496,
       -94.71382002, -94.67915728, -94.64472807, -94.61028893,
       -94.57586149, -94.5413006 , -94.50681276, -94.47232647,
       -94.43783315, -94.40334087, -94.36896439, -94.3344838 ,
       -94.30000226, -94.71379884, -94.67934353, -94.64480863,
       -94.61032437, -94.5758142 , -94.54131858, -94.5068305 ,
       -94.47234373, -94.43785282, -94.40336315, -94.36896752,
       -94.33448358, -94.67931044, -94.64481528, -94.6103127 ,
       -94.57581226, -94.54132421, -94.50686823, -94.47236063,
       -94.43788329, -94.40339566, -94.36896597, -94.33448414,
       -94.6103157 , -94.57582551, -94.5413455 , -94.50686484,
       -94.47238552, -94.43793517, -94.40344988, -94.36896774,
       -94.57586799, -94.54137521, -94.50689182, -94.47241832,
       -94.43794337, -94.40345004, -94.36896517, -94.43793175, -94.40344938])
    lat0 = np.array([ 28.76535526,  28.76708746,  28.77501817,  28.78062341,
        28.78306925,  28.78368291,  28.78405357,  28.78166211,
        28.78210243,  28.7773216 ,  28.77323811,  28.81010872,
        28.81143065,  28.80826667,  28.80367884,  28.79973545,
        28.80031157,  28.80913415,  28.81216795,  28.81724995,
        28.81806324,  28.81905652,  28.81689134,  28.81767776,
        28.81565743,  28.81005498,  28.80735808,  28.80525439,
        28.84540341,  28.84536209,  28.8441454 ,  28.8432604 ,
        28.84498021,  28.84500858,  28.84041811,  28.83543409,
        28.83267942,  28.84392144,  28.84316749,  28.84883649,
        28.8515854 ,  28.85201573,  28.85187791,  28.85077261,
        28.85234814,  28.84776244,  28.84301094,  28.84103623,
        28.83895653,  28.87968787,  28.88097192,  28.88143801,
        28.88028398,  28.87858054,  28.87759731,  28.87899065,
        28.87790163,  28.87165387,  28.86843928,  28.87180208,
        28.87696446,  28.88007685,  28.88400872,  28.88534883,
        28.88571169,  28.88443633,  28.88971935,  28.88451882,
        28.87810314,  28.87425386,  28.87322817,  28.87337358,
        28.87350915,  28.91457525,  28.91536358,  28.91427486,
        28.91469658,  28.91397609,  28.91350777,  28.91200051,
        28.91335948,  28.91095661,  28.90555654,  28.90585106,
        28.91250768,  28.9136333 ,  28.91604794,  28.91683891,
        28.91761613,  28.91760556,  28.92282959,  28.91694153,
        28.91350456,  28.90915139,  28.90835633,  28.90803861,
        28.90878897,  28.90869713,  28.90856713,  28.94995239,
        28.948957  ,  28.94848145,  28.94688026,  28.94898937,
        28.9487707 ,  28.94743282,  28.94719299,  28.94537499,
        28.94232744,  28.94462167,  28.94477507,  28.9493948 ,
        28.95161127,  28.9509314 ,  28.9540917 ,  28.95368218,
        28.94824614,  28.94512331,  28.94378295,  28.94327031,
        28.94317453,  28.94353729,  28.94465235,  28.94466869,
        28.94377094,  28.94363231,  28.98380131,  28.98282002,
        28.98245804,  28.98282346,  28.98569962,  28.98495583,
        28.98377018,  28.98053817,  28.97996426,  28.9791274 ,
        28.9812189 ,  28.98729193,  28.98488803,  28.98536505,
        28.98347972,  28.98314284,  28.97807627,  28.97721417,
        28.97760071,  28.97867585,  28.97937764,  28.97938372,
        28.98050145,  28.98004351,  28.97984485,  29.01903772,
        29.01710296,  29.01603716,  29.0177117 ,  29.02159048,
        29.02168601,  29.02030093,  29.01791038,  29.01628313,
        29.01671817,  29.01695379,  29.01506836,  29.01715138,
        29.01309368,  29.01190686,  29.01212005,  29.0125606 ,
        29.01278356,  29.01396788,  29.01475479,  29.01580351,
        29.01596217,  29.01649367,  29.01561197,  29.05457128,
        29.05200606,  29.05015834,  29.05009941,  29.05569211,
        29.05752021,  29.05555181,  29.05340123,  29.05248823,
        29.05183564,  29.0497549 ,  29.04878452,  29.04821766,
        29.04701212,  29.04736378,  29.04771414,  29.04834269,
        29.04881047,  29.05075742,  29.05181285,  29.0525053 ,
        29.05239873,  29.05188441,  29.08674906,  29.08623775,
        29.08473911,  29.08316682,  29.09122286,  29.09223805,
        29.08911222,  29.08598126,  29.08553778,  29.08482437,
        29.0844318 ,  29.08312266,  29.08156522,  29.08424762,
        29.08400764,  29.08372939,  29.08408783,  29.08503082,
        29.08782164,  29.08828501,  29.08775941,  29.08686735,
        29.11932383,  29.11869667,  29.1200974 ,  29.12473517,
        29.12451615,  29.12164612,  29.12055087,  29.1199654 ,
        29.11966499,  29.11931263,  29.11842251,  29.11853281,
        29.11771036,  29.11753259,  29.11782179,  29.11993195,
        29.12164797,  29.12222157,  29.1217534 ,  29.12120819,
        29.12028068,  29.15294612,  29.15389356,  29.15557098,
        29.15644469,  29.15646748,  29.15550048,  29.1549246 ,
        29.1545665 ,  29.15449671,  29.15455857,  29.15380974,
        29.15268776,  29.15125316,  29.15252149,  29.15265188,
        29.15489481,  29.15578787,  29.15534361,  29.15483907,
        29.15439123,  29.18865454,  29.19040324,  29.19063081,
        29.19045645,  29.19022258,  29.19001559,  29.19009219,
        29.19021768,  29.18894406,  29.18823023,  29.18668596,
        29.1854842 ,  29.18506405,  29.18712063,  29.18821012,
        29.18840919,  29.18854905,  29.18857587,  29.22244242,
        29.22415138,  29.22521163,  29.22503209,  29.22455719,
        29.22443403,  29.2241415 ,  29.22428364,  29.22241002,
        29.22325038,  29.22195349,  29.22105562,  29.22216715,
        29.22263237,  29.22342513,  29.22327686,  29.22285473,
        29.25855078,  29.25948209,  29.25934358,  29.25909356,
        29.25897595,  29.25858898,  29.25800598,  29.25716689,
        29.25854956,  29.25718806,  29.25644393,  29.25712975,
        29.25634803,  29.25707811,  29.25718925,  29.29271887,
        29.29383605,  29.29372184,  29.29352412,  29.29339529,
        29.29318299,  29.29357431,  29.29382093,  29.29303787,
        29.29235026,  29.29172412,  29.29147656,  29.29183109,
        29.3265459 ,  29.32779698,  29.32821755,  29.32835817,
        29.3283959 ,  29.32944633,  29.32940504,  29.32871914,
        29.3281645 ,  29.32751699,  29.32641967,  29.32597603,
        29.3257056 ,  29.3608033 ,  29.36081462,  29.36128103,
        29.36245645,  29.36267607,  29.36296771,  29.36435431,
        29.36406126,  29.36266981,  29.36237506,  29.36211035,
        29.36123055,  29.39600219,  29.39537731,  29.39481935,
        29.39562938,  29.39595476,  29.39667924,  29.39832244,
        29.39839242,  29.39724229,  29.39699688,  29.39761639,
        29.43068238,  29.43005703,  29.42975255,  29.42983812,
        29.42961895,  29.43037096,  29.43089437,  29.43062845,
        29.4643337 ,  29.46461921,  29.46480006,  29.46508309,
        29.46469459,  29.46449752,  29.4640883 ,  29.4996813 ,  29.49898067])

    # # Eliminate points that are outside domain or in masked areas
    # lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)

    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]
    # Change to get positions at the center of the given cell
    lon0, lat0, _ = tracpy.tools.interpolate2d(ia - 0.5, ja - 0.5, grid, 'm_ij2ll')
    N = 1 #lon0.size since there is only one drifter per box in this setup

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
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()
    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = (abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/N

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
        grid = tracpy.inout.readgrid(loc)
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
    ndays = 60
    ff = 1 # This is a backward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    # lon0,lat0 = np.meshgrid(np.linspace(-95.3,-94.3,30), 
    #                         np.linspace(28.6,29.6,30))
    # These are from the first simulation of Matt's, directly (5-23-10-T00)
    lon0 = np.array([-95.0241399 , -94.98965578, -94.95517325, -94.92069183,
       -94.88620883, -94.85172509, -94.81724335, -94.78275982,
       -94.7482779 , -94.71379507, -94.67931158, -95.16206948,
       -95.12758717, -95.09310538, -95.05862145, -95.02413794,
       -94.98965701, -94.95513931, -94.9206318 , -94.88614683,
       -94.85167001, -94.81717435, -94.78272419, -94.7482423 ,
       -94.71379539, -94.67931198, -94.64482931, -94.61034688,
       -95.2655187 , -95.2310361 , -95.19655087, -95.16207005,
       -95.12758782, -95.09305431, -95.05857084, -95.02410113,
       -94.98962946, -94.95511785, -94.92061517, -94.88611981,
       -94.85162388, -94.81713473, -94.78264482, -94.7481572 ,
       -94.71366432, -94.67922921, -94.64482815, -94.61034697,
       -94.57586365, -95.3000017 , -95.2655181 , -95.2310333 ,
       -95.19653484, -95.1620388 , -95.12755681, -95.0930927 ,
       -95.0586025 , -95.02410465, -94.98961709, -94.95512324,
       -94.92061694, -94.88611758, -94.85162163, -94.81713356,
       -94.78263082, -94.74814088, -94.71364736, -94.67914544,
       -94.64466756, -94.61020375, -94.57586434, -94.5413805 ,
       -94.50689871, -95.29999722, -95.26551799, -95.23105674,
       -95.19657054, -95.16205672, -95.12757422, -95.09306828,
       -95.05860208, -95.02411644, -94.9896116 , -94.95511559,
       -94.92060986, -94.8861118 , -94.85162132, -94.81713422,
       -94.78264541, -94.74815372, -94.71370889, -94.67918341,
       -94.64467528, -94.61020068, -94.57572683, -94.54138145,
       -94.50689755, -94.47241587, -94.4379322 , -95.26555061,
       -95.23103205, -95.19653341, -95.16204189, -95.12755044,
       -95.0930626 , -95.0585794 , -95.02409261, -94.98960462,
       -94.95510265, -94.9206037 , -94.88610817, -94.85162115,
       -94.81713221, -94.78264374, -94.74815832, -94.71366663,
       -94.67925681, -94.64469702, -94.61022039, -94.57573748,
       -94.54125526, -94.50679148, -94.47241496, -94.43793323,
       -94.40344925, -94.36896757, -95.19649488, -95.1620103 ,
       -95.12752738, -95.0930437 , -95.05855748, -95.02407384,
       -94.98957927, -94.95508959, -94.92060034, -94.88611059,
       -94.85161738, -94.81712686, -94.78264129, -94.74815407,
       -94.71366744, -94.67918942, -94.64470334, -94.61022465,
       -94.57583904, -94.54130193, -94.50677423, -94.47231023,
       -94.43793233, -94.40345044, -94.36896651, -95.16199311,
       -95.12750233, -95.09302201, -95.05853444, -95.02404783,
       -94.98955697, -94.95507538, -94.9205921 , -94.88610191,
       -94.85161473, -94.81712553, -94.78263937, -94.74815334,
       -94.71366913, -94.67929834, -94.64471377, -94.61023533,
       -94.57575945, -94.54126356, -94.5067751 , -94.47228261,
       -94.43793208, -94.40344962, -94.36896753, -95.1275201 ,
       -95.09300681, -95.05853101, -95.0240449 , -94.989561  ,
       -94.95507202, -94.92058516, -94.88609511, -94.85161244,
       -94.8171244 , -94.78263975, -94.74815675, -94.7136751 ,
       -94.67919383, -94.64470942, -94.61022975, -94.57574723,
       -94.54126392, -94.50677672, -94.47228188, -94.43779522,
       -94.40344873, -94.36896673, -95.09307524, -95.0585704 ,
       -95.02404034, -94.98955968, -94.95506705, -94.9205787 ,
       -94.88609684, -94.85161301, -94.81712623, -94.78264071,
       -94.74816118, -94.71367794, -94.67919103, -94.64470725,
       -94.61022116, -94.57574044, -94.54125494, -94.50676922,
       -94.4722896 , -94.43779573, -94.40330904, -94.36896617,
       -95.0240494 , -94.98955984, -94.95506912, -94.92058844,
       -94.88610242, -94.8516158 , -94.81713094, -94.78264925,
       -94.74816199, -94.71367295, -94.67919152, -94.64470357,
       -94.61021991, -94.57573713, -94.54125231, -94.50676919,
       -94.47227957, -94.43779669, -94.40330223, -94.36896737,
       -94.33448354, -94.98957868, -94.955084  , -94.92059823,
       -94.88611136, -94.85162801, -94.81714058, -94.78265147,
       -94.74816622, -94.71367877, -94.67919459, -94.64470943,
       -94.6102225 , -94.57573871, -94.54125596, -94.50676809,
       -94.47227971, -94.43779329, -94.40330804, -94.36896646,
       -94.33448473, -94.920616  , -94.88612119, -94.85163229,
       -94.81714512, -94.78265761, -94.74816964, -94.71379487,
       -94.67919069, -94.64470996, -94.61022605, -94.5757423 ,
       -94.54125522, -94.50677662, -94.47228688, -94.43779977,
       -94.40331201, -94.36896743, -94.33448378, -94.88620814,
       -94.85165828, -94.81716063, -94.78267083, -94.74824137,
       -94.71368664, -94.67920442, -94.64472177, -94.61022466,
       -94.57574458, -94.54126575, -94.50678072, -94.47229048,
       -94.43780518, -94.40331805, -94.36896525, -94.3344843 ,
       -94.8171858 , -94.78268892, -94.74820967, -94.71371169,
       -94.67922681, -94.64471383, -94.61023142, -94.57575611,
       -94.54126891, -94.50678384, -94.47229335, -94.43781414,
       -94.40331992, -94.36882327, -94.33448285, -94.74824019,
       -94.71373235, -94.67923883, -94.64476681, -94.61023978,
       -94.57575757, -94.54128176, -94.50679793, -94.47231165,
       -94.43781791, -94.40331198, -94.36896724, -94.33448496,
       -94.71382002, -94.67915728, -94.64472807, -94.61028893,
       -94.57586149, -94.5413006 , -94.50681276, -94.47232647,
       -94.43783315, -94.40334087, -94.36896439, -94.3344838 ,
       -94.30000226, -94.71379884, -94.67934353, -94.64480863,
       -94.61032437, -94.5758142 , -94.54131858, -94.5068305 ,
       -94.47234373, -94.43785282, -94.40336315, -94.36896752,
       -94.33448358, -94.67931044, -94.64481528, -94.6103127 ,
       -94.57581226, -94.54132421, -94.50686823, -94.47236063,
       -94.43788329, -94.40339566, -94.36896597, -94.33448414,
       -94.6103157 , -94.57582551, -94.5413455 , -94.50686484,
       -94.47238552, -94.43793517, -94.40344988, -94.36896774,
       -94.57586799, -94.54137521, -94.50689182, -94.47241832,
       -94.43794337, -94.40345004, -94.36896517, -94.43793175, -94.40344938])
    lat0 = np.array([ 28.76535526,  28.76708746,  28.77501817,  28.78062341,
        28.78306925,  28.78368291,  28.78405357,  28.78166211,
        28.78210243,  28.7773216 ,  28.77323811,  28.81010872,
        28.81143065,  28.80826667,  28.80367884,  28.79973545,
        28.80031157,  28.80913415,  28.81216795,  28.81724995,
        28.81806324,  28.81905652,  28.81689134,  28.81767776,
        28.81565743,  28.81005498,  28.80735808,  28.80525439,
        28.84540341,  28.84536209,  28.8441454 ,  28.8432604 ,
        28.84498021,  28.84500858,  28.84041811,  28.83543409,
        28.83267942,  28.84392144,  28.84316749,  28.84883649,
        28.8515854 ,  28.85201573,  28.85187791,  28.85077261,
        28.85234814,  28.84776244,  28.84301094,  28.84103623,
        28.83895653,  28.87968787,  28.88097192,  28.88143801,
        28.88028398,  28.87858054,  28.87759731,  28.87899065,
        28.87790163,  28.87165387,  28.86843928,  28.87180208,
        28.87696446,  28.88007685,  28.88400872,  28.88534883,
        28.88571169,  28.88443633,  28.88971935,  28.88451882,
        28.87810314,  28.87425386,  28.87322817,  28.87337358,
        28.87350915,  28.91457525,  28.91536358,  28.91427486,
        28.91469658,  28.91397609,  28.91350777,  28.91200051,
        28.91335948,  28.91095661,  28.90555654,  28.90585106,
        28.91250768,  28.9136333 ,  28.91604794,  28.91683891,
        28.91761613,  28.91760556,  28.92282959,  28.91694153,
        28.91350456,  28.90915139,  28.90835633,  28.90803861,
        28.90878897,  28.90869713,  28.90856713,  28.94995239,
        28.948957  ,  28.94848145,  28.94688026,  28.94898937,
        28.9487707 ,  28.94743282,  28.94719299,  28.94537499,
        28.94232744,  28.94462167,  28.94477507,  28.9493948 ,
        28.95161127,  28.9509314 ,  28.9540917 ,  28.95368218,
        28.94824614,  28.94512331,  28.94378295,  28.94327031,
        28.94317453,  28.94353729,  28.94465235,  28.94466869,
        28.94377094,  28.94363231,  28.98380131,  28.98282002,
        28.98245804,  28.98282346,  28.98569962,  28.98495583,
        28.98377018,  28.98053817,  28.97996426,  28.9791274 ,
        28.9812189 ,  28.98729193,  28.98488803,  28.98536505,
        28.98347972,  28.98314284,  28.97807627,  28.97721417,
        28.97760071,  28.97867585,  28.97937764,  28.97938372,
        28.98050145,  28.98004351,  28.97984485,  29.01903772,
        29.01710296,  29.01603716,  29.0177117 ,  29.02159048,
        29.02168601,  29.02030093,  29.01791038,  29.01628313,
        29.01671817,  29.01695379,  29.01506836,  29.01715138,
        29.01309368,  29.01190686,  29.01212005,  29.0125606 ,
        29.01278356,  29.01396788,  29.01475479,  29.01580351,
        29.01596217,  29.01649367,  29.01561197,  29.05457128,
        29.05200606,  29.05015834,  29.05009941,  29.05569211,
        29.05752021,  29.05555181,  29.05340123,  29.05248823,
        29.05183564,  29.0497549 ,  29.04878452,  29.04821766,
        29.04701212,  29.04736378,  29.04771414,  29.04834269,
        29.04881047,  29.05075742,  29.05181285,  29.0525053 ,
        29.05239873,  29.05188441,  29.08674906,  29.08623775,
        29.08473911,  29.08316682,  29.09122286,  29.09223805,
        29.08911222,  29.08598126,  29.08553778,  29.08482437,
        29.0844318 ,  29.08312266,  29.08156522,  29.08424762,
        29.08400764,  29.08372939,  29.08408783,  29.08503082,
        29.08782164,  29.08828501,  29.08775941,  29.08686735,
        29.11932383,  29.11869667,  29.1200974 ,  29.12473517,
        29.12451615,  29.12164612,  29.12055087,  29.1199654 ,
        29.11966499,  29.11931263,  29.11842251,  29.11853281,
        29.11771036,  29.11753259,  29.11782179,  29.11993195,
        29.12164797,  29.12222157,  29.1217534 ,  29.12120819,
        29.12028068,  29.15294612,  29.15389356,  29.15557098,
        29.15644469,  29.15646748,  29.15550048,  29.1549246 ,
        29.1545665 ,  29.15449671,  29.15455857,  29.15380974,
        29.15268776,  29.15125316,  29.15252149,  29.15265188,
        29.15489481,  29.15578787,  29.15534361,  29.15483907,
        29.15439123,  29.18865454,  29.19040324,  29.19063081,
        29.19045645,  29.19022258,  29.19001559,  29.19009219,
        29.19021768,  29.18894406,  29.18823023,  29.18668596,
        29.1854842 ,  29.18506405,  29.18712063,  29.18821012,
        29.18840919,  29.18854905,  29.18857587,  29.22244242,
        29.22415138,  29.22521163,  29.22503209,  29.22455719,
        29.22443403,  29.2241415 ,  29.22428364,  29.22241002,
        29.22325038,  29.22195349,  29.22105562,  29.22216715,
        29.22263237,  29.22342513,  29.22327686,  29.22285473,
        29.25855078,  29.25948209,  29.25934358,  29.25909356,
        29.25897595,  29.25858898,  29.25800598,  29.25716689,
        29.25854956,  29.25718806,  29.25644393,  29.25712975,
        29.25634803,  29.25707811,  29.25718925,  29.29271887,
        29.29383605,  29.29372184,  29.29352412,  29.29339529,
        29.29318299,  29.29357431,  29.29382093,  29.29303787,
        29.29235026,  29.29172412,  29.29147656,  29.29183109,
        29.3265459 ,  29.32779698,  29.32821755,  29.32835817,
        29.3283959 ,  29.32944633,  29.32940504,  29.32871914,
        29.3281645 ,  29.32751699,  29.32641967,  29.32597603,
        29.3257056 ,  29.3608033 ,  29.36081462,  29.36128103,
        29.36245645,  29.36267607,  29.36296771,  29.36435431,
        29.36406126,  29.36266981,  29.36237506,  29.36211035,
        29.36123055,  29.39600219,  29.39537731,  29.39481935,
        29.39562938,  29.39595476,  29.39667924,  29.39832244,
        29.39839242,  29.39724229,  29.39699688,  29.39761639,
        29.43068238,  29.43005703,  29.42975255,  29.42983812,
        29.42961895,  29.43037096,  29.43089437,  29.43062845,
        29.4643337 ,  29.46461921,  29.46480006,  29.46508309,
        29.46469459,  29.46449752,  29.4640883 ,  29.4996813 ,  29.49898067])

    # # Eliminate points that are outside domain or in masked areas
    # lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)

    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]
    # Change to get positions at the center of the given cell
    lon0, lat0, _ = tracpy.tools.interpolate2d(ia - 0.5, ja - 0.5, grid, 'm_ij2ll')
    N = 1 #lon0.size since there is only one drifter per box in this setup

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
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()
    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = (abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/N

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
        grid = tracpy.inout.readgrid(loc)
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

    lon0, lat0 = tracpy.tools.seed(lon, lat, dlon=dlon, dlat=dlat, N=Nh)

    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)

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
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    lon0,lat0 = np.meshgrid(np.linspace(-88.81-.4,-88.81+.4,15), 
                            np.linspace(28.1,28.9,15))

    # Eliminate points that are outside domain or in masked areas
    lon0,lat0 = tracpy.tools.check_points(lon0,lat0,grid)

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
        grid = tracpy.inout.readgrid(loc)
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
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]
    # Change to get positions at the center of the given cell
    lon0, lat0, _ = tracpy.tools.interpolate2d(ia - 0.5, ja - 0.5, grid, 'm_ij2ll')

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
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
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
        grid = tracpy.inout.readgrid(loc)
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
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]
    # Change to get positions at the center of the given cell
    lon0, lat0, _ = tracpy.tools.interpolate2d(ia - 0.5, ja - 0.5, grid, 'm_ij2ll')

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
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
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

def gom_dwh_f(date, N, grid=None):
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
    # loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
    # loc = 'http://omgsrv1.meas.ncsu.edu:8080/thredds/dodsC/fmrc/sabgom_tracer_out2/SABGOM_Model_Run_Collection_Tracer_out2_best.ncd'
    # loc = 'http://omgsrv1.meas.ncsu.edu:8080/thredds/dodsC/fmrc/sabgom_tracer_out1/SABGOM_Model_Run_Collection_Tracer_out1_best.ncd'
    loc = 'http://omgsrv1.meas.ncsu.edu:8080/thredds/dodsC/fmrc/sabgom/SABGOM_Forecast_Model_Run_Collection_best.ncd'
    # loc = 'http://omgsrv1.meas.ncsu.edu:8080/thredds/sabgom_catalog.html'
    # loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/fmrc/roms/out/ROMS_Output_Feature_Collection_Aggregation_best.ncd'

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ndays = 60
    ff = 1 # This is a forward-moving simulation

    # Time between outputs
    tseas = 3*3600 # 3 hours between outputs, in seconds, time between model outputs 
    ah = 20.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    # These are the actual blowout coords for DWH (shifted to be in domain)
    lon0 = np.array([-88.36594444444444])
    lat0 = np.array([28.73813888888889])
    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')
    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]
    # Change to get positions at the center of the given cell
    lon0, lat0, _ = tracpy.tools.interpolate2d(ia - 0.5, ja - 0.5, grid, 'm_ij2ll')

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
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
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
        name = 'gom_dwh_f/' + date.isoformat()[0:13] + 'N' + str(N)

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid, dostream, T0, U, V


def allgrid_f(date=None, grid=None):
    '''
    Initialization for seeding drifters at all shelf model grid points to be run
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
    ff = 1 # This is a forward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 0.
    av = 0. # m^2/s

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    # Use the center of all grid cells
    dx = 3; dy = 3;
    lon0 = grid['lonr'][1:-1:dx,1:-1:dy]
    lat0 = grid['latr'][1:-1:dx,1:-1:dy]

    # Eliminate points that are outside domain or in masked areas
    lon0, lat0 = tracpy.tools.check_points(lon0, lat0, grid)

    # Interpolate to get starting positions in grid space
    xstart0, ystart0, _ = tracpy.tools.interpolate2d(lon0, lat0, grid, 'd_ll2ij')

    # Initialize seed locations 
    ia = np.ceil(xstart0).astype(int) #[253]#,525]
    ja = np.ceil(ystart0).astype(int) #[57]#,40]
    # lon0, lat0 already at cell centers
    # # Change to get positions at the center of the given cell
    # lon0, lat0, _ = tracpy.tools.interpolate2d(ia - 0.5, ja - 0.5, grid, 'm_ij2ll')
    N = 1 #lon0.size since there is only one drifter per box in this setup

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
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()
    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = (abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/N

    # Initialize the arrays to save the transports on the grid in the loop.
    # These arrays aggregate volume transport when a drifter enters or exits a grid cell
    # These should start at zero since we don't know which way things will travel yet
    U = np.ma.zeros(grid['xu'].shape,order='F')
    V = np.ma.zeros(grid['xv'].shape,order='F')

    # simulation name, used for saving results into netcdf file
    if date is None:
        name = 'temp' #'5_5_D5_F'
    else:
        name = 'allgrid_f/' + date.isoformat()[0:13] 

    return loc, nsteps, ndays, ff, date, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, name, grid, dostream, T0.data, U, V


def galvcon_b(date, loc, grid=None):
    '''
    Initialization for seeding drifters near Galveston Bay to be run
    backward over a long period of time to look at Bay connectivity.

    Optional inputs for making tests easy to run:
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    # Initialize parameters
    nsteps = 5 # 5 time interpolation steps
    ndays = 30#180 # in days, about 6 months
    ff = -1 # This is a backward-moving simulation

    # Time between outputs
    tseas = 4*3600 # 4 hours between outputs, in seconds, time between model outputs 
    ah = 20.
    av = 0. # m^2/s

    # Can use a subset of the drifters in order to test sensitivity to N
    N = 1000

    if grid is None:
        # if loc is the aggregated thredds server, the grid info is
        # included in the same file
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Initial lon/lat locations for drifters
    # Choosing a point source an each of the two openings to the Bay.
    # Rho indices ia=[257,279] and ja=[152,160], so already at center of cells
    ia0 = (257,279); ja0 = (152,160);
    ia = np.concatenate((np.ones(N/2,np.dtype(int))*ia0[0], np.ones(N/2,np.dtype(int))*ia0[1]))
    ja = np.concatenate((np.ones(N/2,np.dtype(int))*ja0[0], np.ones(N/2,np.dtype(int))*ja0[1]))
    lon0 = grid['lonr'][ia,ja]
    lat0 = grid['latr'][ia,ja]

    # Start at centers in grid space too
    xstart0, ystart0 = ia - 0.5, ja - 0.5

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
    nc, tinds = tracpy.inout.setupROMSfiles(loc, datenum, ff, tout)
    # Get fluxes at first time step in order to find initial drifter volume transport
    uf, vf, dzt, zrt, zwt  = tracpy.inout.readfields(tinds[0],grid,nc,z0,zpar)
    nc.close()
    # Initial total volume transport as a scalar quantity to be conserved, I think
    T0 = ((abs(uf[ia, ja, 0]) + abs(vf[ia, ja, 0]))/(N/2.)).data

    # Initialize the arrays to save the transports on the grid in the loop.
    # These arrays aggregate volume transport when a drifter enters or exits a grid cell
    # These should start at zero since we don't know which way things will travel yet
    U = np.ma.zeros(grid['xu'].shape,order='F')
    V = np.ma.zeros(grid['xv'].shape,order='F')

    return nsteps, ndays, ff, tseas, ah, av, lon0, lat0, \
            z0, zpar, do3d, doturb, grid, dostream, N, T0, U, V
