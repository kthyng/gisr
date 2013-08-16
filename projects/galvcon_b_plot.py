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
ndays = 60
# d drifters for quiver
dd = 6

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

for File in Files:

	# Find tinds for track in full model output set
	track = netCDF.Dataset(File)
	tp = track.variables['tp'][:]
	tstart = find(tp.max() == t)
	# time indices for ndays with 6 outputs per day (4 hour frequency) (5 interpolation steps)
	tinds = np.arange(tstart, tstart-ndays*6*5, -1) #find(tp.max() == t))

	# Plot a wind array from a representative location in the TXLA domain as a wind legend
	# Read in model output. Negative sign since backward in time
	wi = d.variables['sustr'][tinds,jind,iind] # alongshore component of wind stress
	wj = d.variables['svstr'][tinds,jind,iind] # acrossshore component of wind stress
	theta = d.variables['angle'][jind, iind]

	# Rotate model output onto Cartesian axes
	wx = wi*np.cos(theta) - wj*np.sin(theta)
	wy = wi*np.sin(theta) + wj*np.cos(theta)

	# Negative sign since backward in time
	# also smooth
	wx = -smooth(wx, window_len=dd/2)
	wy = -smooth(wy, window_len=dd/2)

	# # Average model output
	# wxm = np.mean(wx,0)
	# wym = np.mean(wy,0)

	lonp = track.variables['lonp'][:,:len(tinds)]
	latp = track.variables['latp'][:,:len(tinds)]
	name = 'galvcon_b/' + str(ndays) + 'days/' + File[17:30]
	tracpy.plotting.tracks(lonp, latp, name, grid)

	# Plot wind arrows
	lonv = np.linspace(-95.2, -88.3, len(wx))
	latv = np.ones(lonv.shape)*25.5
	x0, y0 = grid['basemap'](lonv, latv)
	plt.quiver(x0[::dd], y0[::dd], wx[::dd], wy[::dd], scale=5, color='grey', width=.003, alpha=.8)
	plt.savefig('figures/' + name + 'tracks.png',bbox_inches='tight')
	plt.close()

	track.close()

d.close()

def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    import numpy as np    
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]

# if __name__ == "__main__":
#     run()    