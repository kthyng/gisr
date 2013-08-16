from __future__ import division
import csv
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

def read():

    # Number of drift cards set out
    ndrift = 1770

    # Read in all initial drift card locations
    idnum = np.zeros(ndrift).astype(int)
    startdate = np.zeros(ndrift)*np.nan
    startlat = np.zeros(ndrift)*np.nan
    startlon = np.zeros(ndrift)*np.nan

    with open('projects/driftcards/Deployment_data.txt','rU') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        for i, row in enumerate(reader):
            idnum[i] = row[0]
            startlat[i] = row[2]
            startlon[i] = row[3]

            temp = row[1]
            year = int(temp.split()[0].split('-')[0])
            month = int(temp.split()[0].split('-')[1])
            day = int(temp.split()[0].split('-')[2])
            hour = int(temp.split()[1].split(':')[0])
            if temp.split()[2] == 'PM' and hour is not 12:
                hour = hour + 12
            minutes = int(temp.split()[1].split(':')[1])
            startdate[i] = netCDF.date2num(datetime(year, month, day, hour, minutes), units)

    # ind = idnum.argsort()
    # Find end locations for given driftcards
    enddate = np.zeros(ndrift)*np.nan
    endlat = np.zeros(ndrift)*np.nan
    endlon = np.zeros(ndrift)*np.nan

    with open('projects/driftcards/Recovery_data.txt','rU') as csvfile:
        reader = csv.reader(csvfile,delimiter=',')
        for i, row in enumerate(reader):
            # this is the index that corresponds to the deployment data
            # can use this to keep everything referenced together
            ind = find(idnum==int(row[0]))
            endlat[ind] = row[2]
            endlon[ind] = row[3]

            temp = row[1]
            year = int(temp.split()[0].split('-')[0])
            month = int(temp.split()[0].split('-')[1])
            day = int(temp.split()[0].split('-')[2])
            hour = int(temp.split()[1].split(':')[0])
            if temp.split()[2] == 'PM' and hour is not 12:
                hour = hour + 12
            minutes = int(temp.split()[1].split(':')[1])
            enddate[ind] = netCDF.date2num(datetime(year, month, day, hour, minutes), units)

    # All arrays are referenced to the same ordering of the drifter idnum

    # Find drifters that started within the TX-LA numerical domain

    # For drifters that started within domain, order by start date and time
    # Start drifters at drift card locations, and have a different simulation for
    # each set of drifters by start date/time
    # Run for three months?

    pdb.set_trace()


if __name__ == "__main__":
    read()    