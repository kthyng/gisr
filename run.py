"""
Attempt at having a control-all file for oil modeling paper.
"""

import numpy as np
import sys
import os
import op
import tracmass
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
from matplotlib.pyplot import *
import glob
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import time
import tracpy
import init
from scipy import ndimage

# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

### Sensitivity study (move to separate file later) ###

# Parameters to be rotated through
nsteps = np.array([5, 10, 15])
ahs = np.array([5., 20., 100.])
nlon = np.array([220, 110, 21])
nlat = np.array([196, 98, 20])
nspaces = np.list(['5', '10', '50']) # names for nlon, nlat
doturb = np.array([0, 1, 2, 3])
turbnames = np.list(['None', 'Turb', 'D', 'A'])

# Loop through the tests that we want to investigate
for test in xrange(ntests):

	# loop through drifter spacing names
	for m, nspace in enumerate(nspaces):
		# loop through nsteps
		for nstep in nsteps:
			# loop through turb names
			for n, turbname in enumerate(turbnames):
				# loop through horizontal viscosities
				for p, ah in enumerate(ahs):

					# Add information to name
					name = nspace + '_' + str(nstep) + '_' + turbname + str(ah) '_F'

					# Read in simulation initialization
					if grid in locals(): # don't need to reread grid
						loc, nsteps, ndays, ff, date, tseas, \
							ah, av, lon0, lat0, z0, zpar, do3d, \
							doturb, name = \
							init.test1(loc='local', nsteps=nstep, \
										ah=ah, grid=grid, nlon=nlon[m], \
										nlat=nlat[m], doturb=doturb[n], name=name)
					else: # need to read in grid
						loc, nsteps, ndays, ff, date, tseas, \
							ah, av, lon0, lat0, z0, zpar, do3d, \
							doturb, name = \
							init.test1(loc='local', nsteps=nstep, \
										ah=ah, grid=None, nlon=nlon[m], \
										nlat=nlat[m], doturb=doturb[n], name=name)
	# pdb.set_trace()

	# If the particle trajectories have not been run, run them
	if not os.path.exists('tracks/' + name + '.nc'):
		# TODO: Try to put each simulation on a different core of the current machine, except 1 or 2
		lonp, latp, zp, t, grid = tracpy.run.run(loc, nsteps, ndays, ff, date, \
										tseas, ah, av, lon0, lat0, \
										z0, zpar, do3d, doturb, name)

	else: # if the files already exist, just read them in for plotting
		d = netCDF.Dataset('tracks/' + name + '.nc')
		lonp = d.variables['lonp'][:]
		latp = d.variables['latp'][:]

	# pdb.set_trace()
	ln = lonp.shape[1]
	# Save final locations of drifters for summary origin plots
	lonptemp, \
		latptemp = tracpy.tools.find_final(lonp, latp)

	# Plot tracks
	# pdb.set_trace()
	# tracpy.plotting.tracks(lonp,latp,name,grid=grid)

# 	# Plot final location (by time index) histogram
# 	# tracpy.plotting.hist(lonp,latp,name,grid=grid,which='contour')
# 	# xmin, ymin = grid['basemap'](lonp.min()-.1, latp.min()+.1)
# 	# xmax, ymax = grid['basemap'](lonp.max()+.1, latp.max()+.1)
# 	# tracpy.plotting.hist(lonp,latp,name,grid=grid, \
# 	# 						which='pcolor',bins=(80,80), \
# 	# 						xlims=[xmin, xmax], \
# 	# 						ylims=[ymin, ymax])	

# 	# pdb.set_trace()
# 	# ADD ABILITY TO JUST READ IN TRACKS IF ALREADY DONE

# # pdb.set_trace()
# # Make histogram of all final locations
# tracpy.plotting.hist(lonpsave,latpsave,name,grid=grid,tind='vector', \
# 							which='pcolor',bins=(80,80))


### End of sensitivity study ###

# # Compile tex document with figures in it
# # !pdflatex dwight.tex