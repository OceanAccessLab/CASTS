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

The purpose of this script is to convert the current structre of the CIOOS_ERDDAP
netcdf files into yearly ones containing the following variables; time, latitude,
longitude, temperature, salinity, level, instrument ID, and sounder depth

We take three files from this online source, all originally from Maritime region
AZMP, AZOMP, and ecosystems
These files are updated yearly so we'll be able to download new ones and run the 
same script.

Paths to files and output locations are specific to the source computer, change as needed

'''


##########################################################################################################
##(1) - Isolate and reformat casts according to year

#Import the three relevant datasets
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CIOOS_ERRDAP/data_raw/May_2023/'
path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CIOOS_ERRDAP/data_processed/'


#Define the names of the files
files = np.array(os.listdir(path))
files = files[np.array([i.endswith('.nc') for i in files])]

#Pre-define the years covered and the number of levels
years = np.arange(1998,2022+1)
max_level = 5000	

#Cycle through and import the data
for i in years:
	ds = {}

	#Raw data
	latitude,longitude,depth,time,file_names = {},{},{},{},{}
	project,platform_name,platform_id,cruise_name,cruise_number,station,ID,event_number = {},{},{},{},{},{},{},{}
	temp,saln = {},{}
	#Sorted data
	latitude_1D,longitude_1D,time_1D,file_names_1D = {},{},{},{}
	project_1D,platform_name_1D,platform_id_1D,cruise_name_1D,cruise_number_1D,station_1D,ID_1D,event_number_1D = {},{},{},{},{},{},{},{}
	temp_2D,saln_2D = {},{}

	#Cycle through each files
	for ii in files:

		#Load in the dataset
		ds[ii] = xr.open_dataset(path+ii)

		#Select relevant data from year
		ds[ii] = ds[ii].sel(row = ds[ii]['time.year'].values == i)

		#Determine if data is present, if so, continue
		if ds[ii].row.size == 0:
			del(ds[ii]) #If none is, delete the key
		else:

			#Select data that is within the levels covered
			ds[ii] = ds[ii].sel(row = ds[ii].PRESPR01 < max_level)

			#Import the relevant data
			#Dimensions
			latitude[ii] = ds[ii].latitude.values
			longitude[ii] = ds[ii].longitude.values
			depth[ii] = ds[ii].PRESPR01.values
			time[ii] = ds[ii].time.values
			file_names[ii] = np.tile(ii, ds[ii].time.size)
			project[ii] = ds[ii].project.values
			platform_name[ii] = ds[ii].platform_name.values
			platform_id[ii] = ds[ii].platform_id.values
			cruise_name[ii] = ds[ii].cruise_name.values
			cruise_number[ii] = ds[ii].cruise_number.values
			station[ii] = ds[ii].station.values
			ID[ii] = ds[ii].id.values
			event_number[ii] = ds[ii].event_number.values

			#Variables
			temp[ii] = ds[ii].TEMPP901.values
			saln[ii] = ds[ii].PSALST01.values #Subject to change

			#Close
			ds[ii].close()

			#Now change the structure into the familiar format, time by level
			var_size = np.unique(ID[ii]).size
			temp_2D[ii] = np.full((var_size,max_level),np.nan)
			saln_2D[ii] = np.full((var_size,max_level),np.nan)

			#Also change the latitude, longitude and time into simplified 1D arrays
			latitude_1D[ii] = np.zeros((var_size))
			longitude_1D[ii] = np.zeros((var_size))
			project_1D[ii] = np.zeros((var_size),dtype='<U50')
			platform_name_1D[ii] = np.zeros((var_size),dtype='<U25')
			platform_id_1D[ii] = np.zeros((var_size),dtype='<U25')
			cruise_name_1D[ii] = np.zeros((var_size),dtype='<U100')
			cruise_number_1D[ii] = np.zeros((var_size),dtype='<U25')
			station_1D[ii] = np.zeros((var_size),dtype='<U25')
			ID_1D[ii] = np.zeros((var_size),dtype='<U50')
			event_number_1D[ii] = np.zeros((var_size),dtype='<U25')
			time_1D[ii] = []

			#Cycle through each unique time
			for iii in np.arange(var_size):

				#Find the mean temperature in each level bin 
				temp_2D[ii][iii] = stats.binned_statistic(
					depth[ii][ID[ii] == np.unique(ID[ii])[iii]],
					temp[ii][ID[ii] == np.unique(ID[ii])[iii]],
					'mean',
					bins = np.arange(max_level+1)
					).statistic

				#Find the mean salinity in each level bin 
				saln_2D[ii][iii] = stats.binned_statistic(
					depth[ii][ID[ii] == np.unique(ID[ii])[iii]],
					saln[ii][ID[ii] == np.unique(ID[ii])[iii]],
					'mean',
					bins = np.arange(max_level+1)
					).statistic

				#Fill in the 1D variables, check if constant
				if latitude[ii][ID[ii] == np.unique(ID[ii])[iii]].all():
					latitude_1D[ii][iii] = latitude[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('Latitude not constant for cast.')
					latitude_1D[ii][iii] = np.nan
				if longitude[ii][ID[ii] == np.unique(ID[ii])[iii]].all():
					longitude_1D[ii][iii] = longitude[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('Longitude not constant for cast.')
					longitude_1D[ii][iii] = np.nan

				if np.unique(project[ii][ID[ii] == np.unique(ID[ii])[iii]]).size == 1:
					project_1D[ii][iii] = project[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('Project not constant for cast.')
				if np.unique(platform_name[ii][ID[ii] == np.unique(ID[ii])[iii]]).size == 1:
					platform_name_1D[ii][iii] = platform_name[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('platform_name not constant for cast.')
				if np.unique(platform_id[ii][ID[ii] == np.unique(ID[ii])[iii]]).size == 1:
					platform_id_1D[ii][iii] = platform_id[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('platform_id not constant for cast.')
				if np.unique(cruise_name[ii][ID[ii] == np.unique(ID[ii])[iii]]).size == 1:
					cruise_name_1D[ii][iii] = cruise_name[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('cruise_name not constant for cast.')
				if np.unique(cruise_number[ii][ID[ii] == np.unique(ID[ii])[iii]]).size == 1:
					cruise_number_1D[ii][iii] = cruise_number[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('cruise_number not constant for cast.')
				if np.unique(station[ii][ID[ii] == np.unique(ID[ii])[iii]]).size == 1:
					station_1D[ii][iii] = station[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('station not constant for cast.')
				if np.unique(ID[ii][ID[ii] == np.unique(ID[ii])[iii]]).size == 1:
					ID_1D[ii][iii] = ID[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					place1 = list(np.unique(ID[ii][ID[ii] == np.unique(ID[ii])[iii]]))
					ID_1D[ii][iii] = ', '.join(place1)
				if np.unique(event_number[ii][ID[ii] == np.unique(ID[ii])[iii]]).size == 1:
					event_number_1D[ii][iii] = event_number[ii][ID[ii] == np.unique(ID[ii])[iii]][0]
				else:
					print('event_number not constant for cast.')

				if time[ii][ID[ii] == np.unique(ID[ii])[iii]].all():
					time_1D[ii].append(time[ii][ID[ii] == np.unique(ID[ii])[iii]][0])

			time_1D[ii] = np.array(time_1D[ii])
			file_names_1D[ii] = np.tile(ii, time_1D[ii].size)


	#After isolating and sorting data, concatenate together
	temp_2D = np.vstack([temp_2D[ii] for ii in ds])
	saln_2D = np.vstack([saln_2D[ii] for ii in ds])
	latitude_1D = np.concatenate([latitude_1D[ii] for ii in ds])
	longitude_1D = np.concatenate([longitude_1D[ii] for ii in ds])
	file_names_1D = np.concatenate([file_names_1D[ii] for ii in ds])
	time_1D = np.concatenate([time_1D[ii] for ii in ds])
	project_1D = np.concatenate([project_1D[ii] for ii in ds])
	platform_name_1D = np.concatenate([platform_name_1D[ii] for ii in ds])
	platform_id_1D = np.concatenate([platform_id_1D[ii] for ii in ds])
	cruise_name_1D = np.concatenate([cruise_name_1D[ii] for ii in ds])
	cruise_number_1D = np.concatenate([cruise_number_1D[ii] for ii in ds])
	station_1D = np.concatenate([station_1D[ii] for ii in ds])
	ID_1D = np.concatenate([ID_1D[ii] for ii in ds])
	event_number_1D = np.concatenate([event_number_1D[ii] for ii in ds])

	#Sort the data according to time
	temp_2D = temp_2D[np.argsort(time_1D),:]
	saln_2D = saln_2D[np.argsort(time_1D),:]
	latitude_1D = latitude_1D[np.argsort(time_1D)]
	longitude_1D = longitude_1D[np.argsort(time_1D)]
	file_names_1D = file_names_1D[np.argsort(time_1D)]
	project_1D = project_1D[np.argsort(time_1D)]
	platform_name_1D = platform_name_1D[np.argsort(time_1D)]
	platform_id_1D = platform_id_1D[np.argsort(time_1D)]
	cruise_name_1D = cruise_name_1D[np.argsort(time_1D)]
	cruise_number_1D = cruise_number_1D[np.argsort(time_1D)]
	station_1D = station_1D[np.argsort(time_1D)]
	ID_1D = ID_1D[np.argsort(time_1D)]

	time_1D = time_1D[np.argsort(time_1D)]

	#Save the data as a netcdf file 

	#Set up the .nc file
	nc_out = nc.Dataset(path_output+str(i)+'.nc','w')
	
	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'CIOOS_ERDDAP Yearly Netcdf File' #Temporary title for the .nc file
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
	file_names = nc_out.createVariable('file_names', str, ('time'), zlib=True)
	projects = nc_out.createVariable('project', str, ('time'), zlib=True)
	platform_names = nc_out.createVariable('platform_name', str, ('time'), zlib=True)
	platform_ids = nc_out.createVariable('platform_ids', str, ('time'), zlib=True)
	cruise_names = nc_out.createVariable('cruise_name', str, ('time'), zlib=True)
	cruise_numbers = nc_out.createVariable('cruise_number', str, ('time'), zlib=True)
	stations = nc_out.createVariable('station', str, ('time'), zlib=True)
	IDs = nc_out.createVariable('ID', str, ('time'), zlib=True)

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
	file_names[:] = file_names_1D
	times[:] = time_1D
	projects[:] = project_1D
	platform_names[:] = platform_name_1D
	platform_ids[:] = platform_id_1D
	cruise_names[:] = cruise_name_1D
	cruise_numbers[:] = cruise_number_1D
	stations[:] = station_1D
	IDs[:] = ID_1D

	#Fill 2D structure
	temp[:,:] = temp_2D
	saln[:,:] = saln_2D

	#Convert to time stamps
	time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_1D]
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_level)

	#Save and close the .nc file
	nc_out.close()
	print(str(i)+' done.')











