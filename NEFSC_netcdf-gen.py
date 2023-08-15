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

The purpose of this script is to convert the current structre of the NEFSC
.dat files into yearly ones containing the following variables; time, latitude,
longitude, temperature, salinity, level, instrument ID, and sounder depth

This script is based off the netsc.R script written by Chantelle Layton.

The path inputs and outputs are specific to the source computer, change as necessary

'''


##########################################################################################################
##(1) - Convert files from .dat to netcdf

#Declare the path of the files and determine the names of each
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NEFSC/dat_files/'
path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NEFSC/individual_netcdf_files/'
files = os.listdir(path)

#Define the max level
max_level = 5000

#Cycle through each file
for i in files:

	#Import the relevant data
	file = open(path+i, 'r')
	data = file.readlines()

	#Cut out any short rows
	data_row_size = np.array([len(ii) for ii in data])
	data = list(np.array(data)[data_row_size == 81])

	#Headers are marked by a 1 at the 80th (79th) character location
	flag = [int(ii[80-1]) for ii in data]
	flag = np.array(flag) == 1
	header = np.array(data)[flag]

	#Isolate the relevant info for each header
	country_code = {}
	ship = {}
	latitude_unformat = {}
	longitude_unformat = {}
	year = {}
	month = {}
	day = {}
	time_unformat = {}
	cruise_code = {}
	cast_number = {}
	bottom_depth = {}
	deployment_method = {}
	ctd_serial_number = {}
	station_number = {}
	for ii,value in enumerate(header):
		country_code[ii] = value[0:2] 
		ship[ii] = value[2:4]
		latitude_unformat[ii] = value[4:9]
		longitude_unformat[ii] = value[10:15]
		year[ii] = value[16:20]
		month[ii] = value[20:22]
		day[ii] = value[22:24]
		time_unformat[ii] = value[24:28]
		cruise_code[ii] = value[28:30]
		cast_number[ii] = value[30:33]
		bottom_depth[ii] = value[34:37]
		deployment_method[ii] = value[38:39]
		ctd_serial_number[ii] = value[41:45]
		station_number[ii] = value[69:72] #Original lower index was 63

	#Now isolate for the data in each cast
	temperature = {}
	salinity = {}
	depth = {}
	for ii,value in enumerate(header):

		#First determine the locations of the headers
		header_index = np.where(flag == True)[0]

		#Isolate for the relevant data in the cast, account for last cast case
		if header_index[ii] == header_index[-1]:
			cast = data[header_index[ii]+1:-1]
		else:
			cast = data[header_index[ii]+1:header_index[ii+1]]

		depth[ii] = np.array([int(iii[28:31]) for iii in cast])
		temperature[ii] = np.array([iii[32:36] for iii in cast])
		salinity[ii] = np.array([iii[37:42] for iii in cast])

	#Re-format the data/metadata into a more readable output
	latitude = {}
	for ii in latitude_unformat:
		place1 = latitude_unformat[ii].replace(' ','0')
		latitude[ii] = int(place1[0:2]) +\
		int(place1[2:4])/60 +\
		int(place1[4])*(60/3600) 
	longitude = {}
	for ii in longitude_unformat:
		place1 = longitude_unformat[ii].replace(' ','0')
		longitude[ii] = int(place1[0:2]) +\
		int(place1[2:4])/60 +\
		int(place1[4])*(60/3600) 
		longitude[ii] = longitude[ii]*(-1)
	time = {}
	for ii in time_unformat:
		minutes = np.round(int(time_unformat[ii][2:])*0.01*60).astype(int)
		time[ii] = time_unformat[ii][0:2].replace(' ','0') + ':' + '%.2d' % minutes
	#Create a np.datetime64 record for each cast
	#If two spaces '  ' are present, replace with 19
	time_format = {}
	for ii in time_unformat:
		time_format[ii] = np.datetime64(
			year[ii].replace('  ','19') + '-' +\
			month[ii].replace(' ','0') + '-' +\
			day[ii].replace(' ','0') +\
			'T' + time[ii])

	#Re-format the temperature and salinity into 2D arrays
	temperature_format = {}
	for ii in temperature:
		
		#Turn from strings into floats
		place1 = np.full(temperature[ii].size,np.nan)
		for iii,value in enumerate(temperature[ii]):
			if value.isspace():
				None
			else:
				place1[iii] = float((value[0:2]+'.'+value[2:]).replace(' ',''))
		place1[place1 == 99.999] = np.nan

		#Convert to binned means 
		temperature_format[ii] = stats.binned_statistic(
					depth[ii],
					place1,
					'mean',
					bins = np.arange(max_level+1)
					).statistic

	#Same for salinity
	salinity_format = {}
	for ii in salinity:

		#Turn from strings into floats
		place1 = np.full(salinity[ii].size,np.nan)
		for iii,value in enumerate(salinity[ii]):
			if value.isspace():
				None
			else:
				place1[iii] = float((value[0:2]+'.'+value[2:]).replace(' ',''))
		place1[place1 == 99.999] = np.nan

		#Convert to binned means 
		salinity_format[ii] = stats.binned_statistic(
					depth[ii],
					place1,
					'mean',
					bins = np.arange(max_level+1)
					).statistic

	#After isolating and sorting data, concatenate together
	temperature_format = np.vstack([temperature_format[ii] for ii in temperature])
	salinity_format = np.vstack([salinity_format[ii] for ii in salinity])
	latitude_1D = np.array([latitude[ii] for ii in latitude])
	longitude_1D = np.array([longitude[ii] for ii in longitude])
	time_1D = np.array([time_format[ii] for ii in time_format])
	country_code_1D = np.array([country_code[ii] for ii in country_code])
	ship_1D = np.array([ship[ii] for ii in ship])
	cruise_code_1D = np.array([cruise_code[ii] for ii in cruise_code])
	cast_number_1D = np.array([cast_number[ii] for ii in cast_number])
	bottom_depth_1D = np.array([bottom_depth[ii] for ii in bottom_depth])
	deployment_method_1D = np.array([deployment_method[ii] for ii in deployment_method])
	ctd_serial_number_1D = np.array([ctd_serial_number[ii] for ii in ctd_serial_number])
	station_number_1D = np.array([station_number[ii] for ii in station_number])

	#Sort the data according to time
	temperature_format = temperature_format[np.argsort(time_1D),:]
	salinity_format = salinity_format[np.argsort(time_1D),:]
	latitude_1D = latitude_1D[np.argsort(time_1D)]
	longitude_1D = longitude_1D[np.argsort(time_1D)]
	country_code_1D = country_code_1D[np.argsort(time_1D)]
	ship_1D = ship_1D[np.argsort(time_1D)]
	cruise_code_1D = cruise_code_1D[np.argsort(time_1D)]
	cast_number_1D = cast_number_1D[np.argsort(time_1D)]
	bottom_depth_1D = bottom_depth_1D[np.argsort(time_1D)]
	deployment_method_1D = deployment_method_1D[np.argsort(time_1D)]
	ctd_serial_number_1D = ctd_serial_number_1D[np.argsort(time_1D)]
	station_number_1D = station_number_1D[np.argsort(time_1D)]
	time_1D = time_1D[np.argsort(time_1D)]


	#Finally convert into individual netcdf files
	#Set up the .nc file
	nc_out = nc.Dataset(path_output+i[:-4]+'.nc','w')
	
	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'NEFSC Individual Netcdf File' #Temporary title for the .nc file
	nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
	nc_out.description = 'Output by jonathan.coyne@dfo-mpo.gc.ca'
	nc_out.history = 'Created ' + tt.ctime(tt.time())

	#Create dimensions
	time = nc_out.createDimension('time', None) #use date2 for this
	level = nc_out.createDimension('level', max_level) 

	#Create coordinate variables
	times = nc_out.createVariable('time', np.float64, ('time',))
	levels = nc_out.createVariable('level', np.int32, ('level',))

	#Create 1D variables
	latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
	longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
	country_code = nc_out.createVariable('country_code', str, ('time'), zlib=True)
	ship = nc_out.createVariable('ship_code', str, ('time'), zlib=True)
	cruise_code = nc_out.createVariable('cruise_code', str, ('time'), zlib=True)
	cast_number = nc_out.createVariable('cast_number', str, ('time'), zlib=True)
	bottom_depth = nc_out.createVariable('station_bottom_depth', str, ('time'), zlib=True)
	deployment_method = nc_out.createVariable('deployment_method', str, ('time'), zlib=True)
	ctd_serial_number = nc_out.createVariable('CTD_serial_number', str, ('time'), zlib=True)
	station_number = nc_out.createVariable('station_number', str, ('time'), zlib=True)

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
	latitudes[:] = latitude_1D
	longitudes[:] = longitude_1D
	times[:] = time_1D
	country_code[:] = country_code_1D
	ship[:] = ship_1D
	cruise_code[:] = cruise_code_1D
	cast_number[:] = cast_number_1D
	bottom_depth[:] = bottom_depth_1D
	deployment_method[:] = deployment_method_1D
	ctd_serial_number[:] = ctd_serial_number_1D
	station_number[:] = station_number_1D

	#Fill 2D structure
	temp[:,:] = temperature_format
	saln[:,:] = salinity_format

	#Convert to time stamps
	time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_1D]
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_level)

	#Save and close the .nc file
	nc_out.close()
	print(str(i)+' done.')


##########################################################################################################
##(2) - Bring together individual files according to year

#Now gather up all like dat files (in netcdf) and compile into yearly
files = np.array(os.listdir(path_output))
path_yearly_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NEFSC/yearly_netcdf_files/'


#Determine the years of each, then format
file_years = np.zeros(files.size)
for i,value in enumerate(files):
	if value.endswith('.nc'):
		ds = xr.open_dataset(path_output+value)
		if np.unique(ds['time.year']).size == 1:
			file_years[i] = ds['time.year'].values[0]
file_years = file_years.astype(int)


#Cycle through by year
for year in np.unique(file_years):

	#Make sure it's not the '0' year
	if year != 0:

		#Run through each of the files for that year
		ds = {}
		for i in np.where(file_years == year)[0]:

			#Open the files
			ds[i] = xr.open_dataset(path_output+files[i])
			ds[i]['file_names'] = (('time'), np.tile(files[i][:-3]+'.dat',ds[i].time.size))

		#Concatenate together
		ds = xr.concat(list(ds.values()),dim='time')

		#Sort according to time
		ds = ds.sortby('time')

		#Check to see that all time steps are in the year
		if (ds['time.year'].values == year).all():
			print('All steps within '+str(year))
		else:
			print('Mis-matched years, go back and check '+str(year))

			#If mis-matched,isolate for the proper year
			ds = ds.sel(time=ds['time.year'] == year)
			print('Now matching!')

		#Save the file, ensure proper time units
		ds.to_netcdf(path_yearly_output + str(year) + '.nc',encoding={'time':{'units': "seconds since 1900-01-01 00:00:00"}})
		ds.close()



