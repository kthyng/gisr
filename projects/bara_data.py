import csv
import pdb
import numpy as np
from datetime import datetime, timedelta
import netCDF4 as netCDF

'''
Times are in UTC
'''

def retrieve():
	units = 'seconds since 1970-01-01'
	lon = np.zeros((10090,1))*np.nan
	lat = np.zeros((10090,1))*np.nan
	t = []
	month = np.zeros((10090,1))*np.nan
	day = np.zeros((10090,1))*np.nan
	year = np.zeros((10090,1))*np.nan
	hour = np.zeros((10090,1))*np.nan
	minute = np.zeros((10090,1))*np.nan
	date = np.zeros((10090,1))
	j = 0
	with open('Union of All Datasets (Standardized)- Barataria Bay.csv','rb') as csvfile:
		reader = csv.reader(csvfile,delimiter=',')
		for i, row in enumerate(reader):
			if i != 0: #skip the header
				# skip all rows with the following list of descriptors
				# since they aren't oil
				if set(['Biological Index', \
					'Fluorescence Index', \
					'Total Fluorescence', \
					'UV/Visible Light Absorbance at 254 nm', \
					'Temperature',\
					'Dissolved Organic Carbon Concentration', \
					'Density']).isdisjoint(row):
					# some rows are missing info
					if len(row) > 4:
						# only want surface data
						if row[9] == '0':
							# pdb.set_trace()
							lon[j] = row[8]
							lat[j] = row[7]
							t.append(row[10])
							month[j] = int(t[j].split(' ')[0].split('/')[0])
							day[j] = int(t[j].split(' ')[0].split('/')[1])
							year[j] = int(t[j].split(' ')[0].split('/')[2])
							hour[j] = int(t[j].split(' ')[1].split(':')[0])
							minute[j] = int(t[j].split(' ')[1].split(':')[1])
							if t[j].split(' ')[2] == 'PM' and hour[j] != 12:
								hour[j] = hour[j] + 12
							date[j] = netCDF.date2num(datetime(year[j], month[j], day[j], \
												hour[j], minute[j], 0),units)
							j = j + 1

	# eliminate extra nan entries
	ind = ~(np.isnan(lon))
	lon = lon[ind]
	lat = lat[ind]
	month = month[ind]
	day = day[ind]
	year = year[ind]
	hour = hour[ind]
	minute = minute[ind]
	date = date[ind]

	# Sort entries since are a little mixed up
	I = date.argsort() # sort by date
	date = date[I]
	lon = lon[I]
	lat = lat[I]
	month = month[I]
	day = day[I]
	year = year[I]
	hour = hour[I]
	minute = minute[I]

	# Eliminate duplicates
	j = 0
	for i in xrange(len(lon)-1):
		# pdb.set_trace()
		if(lon[i]==lon[i+1] and lat[i]==lat[i+1] and date[i]==date[i+1]):
			lon[i] = np.nan
			lat[i] = np.nan
			t.pop(j)
			month[i] = np.nan
			year[i] = np.nan
			day[i] = np.nan
			hour[i] = np.nan
			minute[i] = np.nan
			date[i] = np.nan
		else:
			j = j + 1

	# eliminate extra nan entries
	ind = ~(np.isnan(lon))
	lon = lon[ind]
	lat = lat[ind]
	month = month[ind]
	day = day[ind]
	year = year[ind]
	hour = hour[ind]
	minute = minute[ind]
	date = date[ind]

	return lon, lat, date