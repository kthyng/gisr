"""
Attempt at having a control-all file for oil modeling paper.
"""

import matplotlib
matplotlib.use("Agg") # set matplotlib to use the backend that does not require a windowing system
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
import projects

# Units for time conversion with netCDF.num2date and .date2num
units = 'seconds since 1970-01-01'

# What tests to run, 1 to run or 0 to skip:
do_sensitivity = 0
do_outer_f = 0
do_galv_b = 0
do_bara_b = 0
do_dwh_f = 1
do_compile = 0


# Run sensitivity tests
if do_sensitivity:
    projects.sensitivity.run()

# Run drifters foward from near the outer numerical boundary
if do_outer_f:
    projects.outer_f.run()

# Run drifters backward from near Galveston Bay
if do_galv_b:
    projects.galv_b.run()

# Run drifters backward from near Barataria Bay
if do_bara_b:
    projects.bara_b.run()

# Run drifters forward from near the Deepwater Horizon Blowout
if do_dwh_f:
    projects.dwh_f.run()



# Compile tex document with figures in it. 
# Run twice to get references correct.
if do_compile:
    os.system("/usr/texbin/pdflatex gisr.tex")
    os.system("/usr/texbin/pdflatex gisr.tex")
