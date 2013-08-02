'''
Functionality to plot drifter transport given the appropriate information
from tracmass.
'''


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
from mpl_toolkits.basemap import Basemap
import op

units = 'seconds since 1970-01-01'

def load(name,fmod=None):
	'''

	Inputs:
		name 	Name of project
		fmod	File modifier: a way to choose a subset of the file in 
				the project directory instead of all. Should be a string and
				can include asterisks as wildcards.

	Outputs:
		U, V 	Transport of drifter volume in x and y directions over all
				used simulation files
		lon0 	Initial lon location for drifters
		lat0 	Initial lat location for drifters
		T0 		Overall
	'''

	# Which files to read in.
	if fmod is None:
		Files = glob.glob('tracks/' + name + '/*.nc')
	else:
		Files = glob.glob('tracks/' + name + '/' + fmod + '.nc')

	Files.sort()

	# Load in U and V volume transports of drifters and add together for
	# all files
	for i, File in enumerate(Files):
	    d = netCDF.Dataset(File)
	    # get lon0, lat0, and initial volumes from init file
	    _, _, _, _, _, _, _, _, \
		    lon0, lat0, _, _, _, _, \
		    _, _, _, T0temp, \
		    _, _ = init.dwh_stream_f(date, 100)

	    if i == 0: # initialize U and V transports from first file
	        U = d.variables['U'][:]
	        V = d.variables['V'][:]
	        T0 = T0temp
	    else: # add in transports from subsequent simulations
	        U = U + d.variables['U'][:]
	        V = V + d.variables['V'][:]
			T0 = T0 + T0temp
	    date = netCDF.num2date(d.variables['tp'][0],units)
	    d.close()

	# # Add initial drifter location (all drifters start at the same location)
	# lon0 = d.variables['lonp'][0,0]
	# lat0 = d.variables['latp'][0,0]

	# old streamline code
	# # Calculate lagrangian barotropic stream function
	# # Doos. psi_i and psi_j should be similar
	# psi_i = np.cumsum(-V, axis=0)[:-1,:]
	# psi_j = np.fliplr(np.cumsum(U[:,::-1], axis=1))[:,:-1]
	# # internet
	# psi_i2 = np.cumsum(V, axis=0)
	# psi_j2 = np.cumsum(U, axis=1)
	# # psi_j2 = np.fliplr(np.cumsum(U[:,::-1], axis=1))
	# psi = psi_j2[:,:-1] - psi_i2[:-1,:]

	return U, V, lon0, lat0, T0

def plot(name, U, V, lon0, lat0, T0, extraname=None):
	'''
	Make plot of zoomed-in area near DWH spill of transport of drifters over 
	time.

	Inputs:
		name
		U
		V
		lon0
		lat0
		T0
	'''

	# Smaller basemap parameters.
	llcrnrlon=-93.5; llcrnrlat=27.2; urcrnrlat=30.7
	loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
	grid = tracpy.inout.readgrid(loc, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, 
	                                urcrnrlat=urcrnrlat)

	fig = plt.figure(figsize=(16.0375,   9.9125))
	tracpy.plotting.background(grid=grid)
	plt.contourf(grid['xpsi'], grid['ypsi'], np.sqrt(op.resize(U,1)**2+op.resize(V,0)**2)/T0, 
	        cmap='gray_r', levels=np.linspace(0,18000,10), extend='max')

	# Inlaid colorbar
	cax = fig.add_axes([0.42, 0.2, 0.27, 0.02])
	cb = colorbar(cax=cax,orientation='horizontal')
	cb.set_label('Normalized drifter transport')
	plt.colorbar()
	plt.title('Deepwater Horizon Spill Transport')

	# Add initial drifter location (all drifters start at the same location)
	x0, y0 = grid['basemap'](lon0, lat0)
	plt.plot(x0, y0, 'go', markersize=10)
	if extraname is None:
		plt.savefig('figures/' + name + '/transport', bbox_inches='tight')
	else:
		plt.savefig('figures/' + name + '/transport' + extraname, bbox_inches='tight')

	# old streamline code
	# # Make plot
	# plt.figure(figsize=(16.0375,   9.9125))
	# contour(grid['xpsi'], grid['ypsi'], abs(psi_i)/abs(psi_i).max(), 20,colors='k',zorder=0)
	# tracpy.plotting.background(grid=grid)
	# # contourf(grid['xr'], grid['yr'], psi, cmap='RdBu_r',
	# #         levels=np.linspace(-600,600,10), 
	# #         extend='both')
	# # contourf(grid['xpsi'], grid['ypsi'], psi_i, cmap='RdBu_r',
	# #         levels=np.linspace(-150000,150000,12), 
	# #         extend='both')
	# # contour(grid['xu'], grid['yu'], psi_j, 0)
	# # colorbar()
	# title('psi_i')
	# # Add initial drifter location (all drifters start at the same location)
	# lon0 = d.variables['lonp'][0,0]
	# lat0 = d.variables['latp'][0,0]
	# x0, y0 = grid['basemap'](lon0, lat0)
	# plot(x0, y0, 'go', markersize=10)
	# # plt.savefig('figures/dwh_stream_f/stream',bbox_inches='tight')


def run():
# def run(name,fmod=None, extraname=None):
	''' Controls which project to run this for'''

	name = 'dwh_stream_f'
	U, V, lon0, lat0, T0 = load(name)
	plot(name, U, V, lon0, lat0, T0)
	# # Load in information
	# if fmod is None:
	# 	U, V, lon0, lat0 = load(name)
	# else:
	# 	U, V, lon0, lat0 = load(name,fmod=fmod)

	# # Plot information
	# if extraname is None:
	# 	plot(name, U, V, lon0, lat0)
	# else:
	# 	plot(name, U, V, lon0, lat0, extraname=extraname)

# def run_dwh_stream_f():
# 	run('dwh_stream_f')

# def run_bara_stream_b():
# 	run('bara_stream_b')

if __name__ == "__main__":
    run()