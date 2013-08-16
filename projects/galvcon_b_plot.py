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
from matplotlib.mlab import find

Files = glob.glob('tracks/galvcon_b/*.nc')
Files.sort()

# number of days to look at
ndays = 10

# Choose a particular wind location
iind= 330; jind = 90;

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
	
d = netCDF.Dataset(loc)
# angle of the grid at the wind point stays the same in time
theta = d.variables['angle'][jind, iind]
# # Interpolate theta to be on psi grid
# theta = op.resize(op.resize(theta,0),1)

# model times
t = d.variables['ocean_time'][:]

grid = tracpy.inout.readgrid(loc)

for File in Files[0:3]:

	# Find tinds for track in full model output set
	track = netCDF.Dataset(File)
	tp = track.variables['tp'][:]
	tstart = find(tp.min() == t)[-1]
	# time indices for ndays with 6 outputs per day (4 hour frequency)
	tinds = np.arange(tstart, tstart-ndays*6, -1) #find(tp.max() == t))

	# Plot a wind array from a representative location in the TXLA domain as a wind legend
	# Read in model output
	wi = d.variables['sustr'][tinds,jind,iind] # alongshore component of wind stress
	wj = d.variables['svstr'][tinds,jind,iind] # acrossshore component of wind stress
	theta = d.variables['angle'][jind, iind]

	# Rotate model output onto Cartesian axes
	wx = wi*np.cos(theta) - wj*np.sin(theta)
	wy = wi*np.sin(theta) + wj*np.cos(theta)

	# # Average model output
	# wxm = np.mean(wx,0)
	# wym = np.mean(wy,0)

	lonp = track.variables['lonp'][:]
	latp = track.variables['latp'][:]
	name = 'galvcon_b/' + File[17:30]
	tracpy.plotting.tracks(lonp, latp, name, grid)

	# Plot wind arrows
	lonv = np.linspace(-95.2, -88.3, len(wx))
	latv = np.ones(lonv.shape)*25.5
	x0, y0 = grid['basemap'](lonv, latv)
	plt.quiver(x0, y0, wx, wy, scale=7, color='grey', width=.003, alpha=.8)
	plt.savefig('figures/' + name + 'tracks.png',bbox_inches='tight')

	track.close()

d.close()

# if __name__ == "__main__":
#     run()    