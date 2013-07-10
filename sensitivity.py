'''
Script to run all steps in the senstivity study part of the GISR work
This runs through a series of parameters that could be important
when running a numerical drifter study, saves the tracks, and plots
the results as tracks and as a histogram.
'''

import numpy as np
import sys
import os
import op
import tracmass
import netCDF4 as netCDF
from mpl_toolkits.basemap import Basemap
import pdb
from matplotlib import delaunay
import matplotlib.pyplot as plt
import glob
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import time
import tracpy
import init
from scipy import ndimage

def run():
	# Make sure necessary directories exist
	if not os.path.exists('tracks'):
		os.makedirs('tracks')
	if not os.path.exists('tracks/sensitivity'):
		os.makedirs('tracks/sensitivity')
	if not os.path.exists('figures'):
		os.makedirs('figures')
	if not os.path.exists('figures/sensitivity'):
		os.makedirs('figures/sensitivity')

	# Parameters to be rotated through
	nsteps = np.array([5, 10, 15])
	ahs = np.array([5, 20])
	nlon = np.array([220, 110, 21])
	nlat = np.array([196, 98, 20])
	nspaces = list(['5', '10', '50']) # names for nlon, nlat
	doturbs = np.array([0, 1, 2, 3])
	turbnames = list(['None', 'Turb', 'D', 'A'])

	# loop through drifter spacing names
	for m, nspace in enumerate(nspaces):
		# loop through nsteps
		for nstep in nsteps:
			# loop through turb names
			for n, turbname in enumerate(turbnames):

				if turbname == 'None': # don't need horizontal viscosity (the input one is not used)

					# Add information to name
					name = 'sensitivity/' + nspace + '_' + str(nstep) + '_' + turbname + '_F'

					# pdb.set_trace()
					# Read in simulation initialization
					# hgrid doesn't have any vertical grid info
					if 'hgrid' in locals(): # don't need to reread grid
						loc, nstep, ndays, ff, date, tseas, \
							ah, av, lon0, lat0, z0, zpar, do3d, \
							doturb, name, hgrid = \
							init.sensitivity(loc='thredds', nsteps=nstep, \
										ah=0, grid=hgrid, nlon=nlon[m], \
										nlat=nlat[m], doturb=doturbs[n], name=name)
					else: # need to read in grid
						loc, nstep, ndays, ff, date, tseas, \
							ah, av, lon0, lat0, z0, zpar, do3d, \
							doturb, name, hgrid = \
							init.sensitivity(loc='thredds', nsteps=nstep, \
										ah=0, grid=None, nlon=nlon[m], \
										nlat=nlat[m], doturb=doturbs[n], name=name)

					# If the particle trajectories have not been run, run them
					if not os.path.exists('tracks/' + name + '.nc'):
						# TODO: Try to put each simulation on a different core of the current machine, except 1 or 2
						if 'grid' in locals(): # don't need to reread grid
							lonp, latp, zp, t, grid = tracpy.run.run(loc, nstep, ndays, ff, date, \
															tseas, ah, av, lon0, lat0, \
															z0, zpar, do3d, doturb, name, grid)
						else:
							lonp, latp, zp, t, grid = tracpy.run.run(loc, nstep, ndays, ff, date, \
															tseas, ah, av, lon0, lat0, \
															z0, zpar, do3d, doturb, name)

					else: # if the files already exist, just read them in for plotting
						d = netCDF.Dataset('tracks/' + name + '.nc')
						lonp = d.variables['lonp'][:]
						latp = d.variables['latp'][:]


					# If the particle trajectories have not been plotted, plot them
					if not os.path.exists('figures/' + name + 'tracks.png'):
						# Plot tracks
						tracpy.plotting.tracks(lonp,latp,name,grid=hgrid)
					if not os.path.exists('figures/' + name + 'histpcolor.png'):
						tracpy.plotting.hist(lonp,latp,name,grid=hgrid, \
											which='pcolor',bins=(40,40))

				else: # need horizontal viscosity
					# loop through horizontal viscosities
					for p, ah in enumerate(ahs):

						# Add information to name
						name = 'sensitivity/' + nspace + '_' + str(nstep) + '_' + turbname + '_F'

						# pdb.set_trace()
						# Read in simulation initialization
						# hgrid doesn't have any vertical grid info
						if 'hgrid' in locals(): # don't need to reread grid
							loc, nstep, ndays, ff, date, tseas, \
								ah, av, lon0, lat0, z0, zpar, do3d, \
								doturb, name, hgrid = \
								init.sensitivity(loc='thredds', nsteps=nstep, \
											ah=ah, grid=hgrid, nlon=nlon[m], \
											nlat=nlat[m], doturb=doturbs[n], name=name)
						else: # need to read in grid
							loc, nstep, ndays, ff, date, tseas, \
								ah, av, lon0, lat0, z0, zpar, do3d, \
								doturb, name, hgrid = \
								init.sensitivity(loc='thredds', nsteps=nstep, \
											ah=ah, grid=None, nlon=nlon[m], \
											nlat=nlat[m], doturb=doturbs[n], name=name)


						# If the particle trajectories have not been run, run them
						if not os.path.exists('tracks/' + name + '.nc'):
							# TODO: Try to put each simulation on a different core of the current machine, except 1 or 2
							if 'grid' in locals(): # don't need to reread grid
								lonp, latp, zp, t, grid = tracpy.run.run(loc, nstep, ndays, ff, date, \
																tseas, ah, av, lon0, lat0, \
																z0, zpar, do3d, doturb, name, grid)
							else:
								lonp, latp, zp, t, grid = tracpy.run.run(loc, nstep, ndays, ff, date, \
																tseas, ah, av, lon0, lat0, \
																z0, zpar, do3d, doturb, name)

						else: # if the files already exist, just read them in for plotting
							d = netCDF.Dataset('tracks/' + name + '.nc')
							lonp = d.variables['lonp'][:]
							latp = d.variables['latp'][:]


						# If the particle trajectories have not been plotted, plot them
						if not os.path.exists('figures/' + name + 'tracks.png'):
							# Plot tracks
							tracpy.plotting.tracks(lonp,latp,name,grid=hgrid)
						if not os.path.exists('figures/' + name + 'histpcolor.png'):
							tracpy.plotting.hist(lonp,latp,name,grid=hgrid, \
												which='pcolor',bins=(40,40))

