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

units = 'seconds since 1970-01-01'

Files = glob.glob('tracks/galvcon_b/*.nc')
Files.sort()

# number of days to look at
ndays = 10
# d drifters for quiver
dd = 1

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
dates = netCDF.num2date(t,units)

# Relative times for outputs
trel = (t-t[0])/(3600.*24) # output times in days

grid = tracpy.inout.readgrid(loc)

for File in Files:

	# Find tinds for track in full model output set
	track = netCDF.Dataset(File)
	tp = track.variables['tp'][:]
	tstart = find(tp.max() == t)
	# time indices for ndays with 6 outputs per day (4 hour frequency)
	tinds_model = np.arange(tstart, tstart-ndays*6, -1) #find(tp.max() == t))
	# time indices for ndays with 6 outputs per day (4 hour frequency) (5 interpolation steps)
	tinds_tracks = np.arange(tstart, tstart-ndays*6*5, -1) #find(tp.max() == t))

	# Plot a wind array from a representative location in the TXLA domain as a wind legend
	# Read in model output. Negative sign since backward in time
	wi = d.variables['sustr'][tinds_model,jind,iind] # alongshore component of wind stress
	wj = d.variables['svstr'][tinds_model,jind,iind] # acrossshore component of wind stress
	theta = d.variables['angle'][jind, iind]

	# Rotate model output onto Cartesian axes
	wx = wi*np.cos(theta) - wj*np.sin(theta)
	wy = wi*np.sin(theta) + wj*np.cos(theta)

	# Negative sign since backward in time
	# also smooth
	wx = -wx
	wy = -wy

	# # Average model output
	# wxm = np.mean(wx,0)
	# wym = np.mean(wy,0)

	lonp = track.variables['lonp'][:,:len(tinds_tracks)]
	latp = track.variables['latp'][:,:len(tinds_tracks)]
	name = 'galvcon_b/' + str(ndays) + 'days/' + File[17:30]
	tracpy.plotting.tracks(lonp, latp, name, grid)

	# Plot wind arrows
	lonv = np.linspace(-95.2, -88.3, len(wx))
	latv = np.ones(lonv.shape)*25.5
	x0, y0 = grid['basemap'](lonv, latv)
	# Plot start and end indicators
	plt.plot(x0[0], y0[0], 'og', markersize=16, alpha=0.5)
	plt.plot(x0[-1], y0[-1], 'or', markersize=16, alpha=0.5)
	# Plot a black line every day on the wind plot
	pdb.set_trace()
	ind = (np.mod(trel[tinds_model],1) == 0.)
	plt.plot(x0[ind], y0[ind], 'k|', markersize=10, alpha=0.5)
	# Plot arrows
	plt.quiver(x0[::dd], y0[::dd], wx[::dd], wy[::dd], scale=5, color='grey', width=.003, alpha=.8)
	# Plot date below wind
	plt.text(x0[ind], y0[ind]-2000, dates[tinds_model][ind])
	plt.savefig('figures/' + name + 'tracks.png',bbox_inches='tight')
	plt.close()

	track.close()

d.close()

# if __name__ == "__main__":
#     run()    