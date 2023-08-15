import xarray as xr
import netCDF4 as nc
import numpy as np
import time as tt
import os
import pfile_tools as p
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt

'''

The purpose of this script is to convert the current structure of the ASCII
netcdf files into yearly ones containing the following variables; time, latitude,
longitude, temperature, salinity, level, instrument ID, and sounder depth

Data was provided by Jean-Luc and Peter Galbraith at MLI, Quebec Region
Updates are provided each year
Casts are included from multiple monitoring programs (AZMP, Ecosystems, etc.)

Casts were provided in one ragged netcdf file. 

'''


##########################################################################################################
##(1) - Isolate and reformat casts according to year

#Import the relevant dataset
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/ASCII/data_raw/'
path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/ASCII/data_processed/'

#Import the raw data
ds = xr.open_dataset(path+'climate_ae.nc')

#Determine the year range covered
years = np.unique(ds['time.year'].values)
max_depth = 5000

#Determine the cast start location for each 
cast_start = np.concatenate((np.array([0]),np.cumsum(ds.cast_size.values)[:-1]))

#Set up empty temperature and salinity arrays
temp_2D = np.full((ds.time.size,max_depth),np.nan)
saln_2D = np.full((ds.time.size,max_depth),np.nan)

#Load in the temperature, salinity, pressure
ds_temp = ds.temperature.values
ds_saln = ds.salinity.values
ds_pres = ds.pressure.values

#Cycle through each of the casts
for i,value in enumerate(cast_start):

	#Slice the temperature, pressure, and salinity
	if i+1 != cast_start.size:
		temp_slice = ds_temp[value:cast_start[i+1]]
		saln_slice = ds_saln[value:cast_start[i+1]]
		pres_slice = ds_pres[value:cast_start[i+1]]
	else:
		ending_index = value + ds.temperature_row_size[-1].values
		temp_slice = ds_temp[value:ending_index]
		saln_slice = ds_saln[value:ending_index]
		pres_slice = ds_pres[value:ending_index]	

	#If the cast is empty, skip
	if pres_slice.size == 0:
		continue

	#Populate the 2D arrays by binning first
	temp_slice = stats.binned_statistic(
		pres_slice,
		temp_slice,
		'mean',
		bins = np.arange(max_depth+1)
		).statistic
	saln_slice = stats.binned_statistic(
		pres_slice,
		saln_slice,
		'mean',
		bins = np.arange(max_depth+1)
		).statistic

	#Populate the 2D arrays
	temp_2D[i,:] = temp_slice
	saln_2D[i,:] = saln_slice 

	print(str(i)+' out of '+str(cast_start.size)+' done.')

#Remove the -99 values, cutoffs 
temp_2D[temp_2D <= -2] = np.nan
temp_2D[temp_2D > 35] = np.nan
saln_2D[saln_2D <= 0] = np.nan
saln_2D[saln_2D > 45] = np.nan

#Record the latitude, longitude, time
years = ds['time.year'].values
ds_time = ds.time.values
ds_latitude = ds.lat.values
ds_longitude = ds.lon.values

#Cycle through each of the years now
for year in np.unique(years):

	#Create a slice of the data for a specific year
	time_slice = ds_time[years == year]
	lat_slice = ds_latitude[years == year]
	lon_slice = ds_longitude[years == year]
	temp_slice = temp_2D[years == year,:]
	saln_slice = saln_2D[years == year,:]

	#Create a filter
	filt = np.argsort(time_slice)

	#Save the data as a netcdf file 
	#Set up the .nc file
	nc_out = nc.Dataset(path_output+str(year)+'.nc','w')
	
	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'MLI Yearly Netcdf File' #Temporary title for the .nc file
	nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
	nc_out.description = 'Output by jonathan.coyne@dfo-mpo.gc.ca'
	nc_out.history = 'Created ' + tt.ctime(tt.time())

	#Create dimensions
	time = nc_out.createDimension('time', None) #use date2 for this
	level = nc_out.createDimension('level', max_depth) 

	#Create coordinate variables
	times = nc_out.createVariable('time', np.float64, ('time',))
	levels = nc_out.createVariable('level', np.int32, ('level',))

	#Create 1D variables
	latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
	longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
	file_names = nc_out.createVariable('file_names',str,('time'),zlib=True)

	#Create 2D variables
	temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	saln = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

	#Variable Attributes
	latitudes.units = 'degree_north'
	longitudes.units = 'degree_east'
	times.units = 'seconds since 1900-01-01 00:00:00'
	times.calendar = 'gregorian'
	levels.units = 'dbar'
	levels.standard_name = "pressure"
	#levels.valid_range = np.array((0.0, 5000.0))
	levels.valid_min = 0
	temp.units = 'Celsius'
	temp.long_name = "Water Temperature" # (may be use to label plots)
	temp.standard_name = "sea_water_temperature"
	saln.long_name = "Practical Salinity"
	saln.standard_name = "sea_water_salinity"
	saln.units = "1"
	saln.valid_min = 0

	#Fill in the 1D variables 
	latitudes[:] = lat_slice[filt]
	longitudes[:] = lon_slice[filt]
	times[:] = time_slice[filt]
	file_names[:] = np.tile('climate_ae.nc', time_slice.size)

	#Fill 2D structure
	temp[:,:] = temp_slice[filt]
	saln[:,:] = saln_slice[filt]

	#Convert to time stamps
	time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_slice[filt]]
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_depth)

	#Save and close the .nc file
	nc_out.close()

	print(str(year)+' done.')








