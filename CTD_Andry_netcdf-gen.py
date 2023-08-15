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

The purpose of this script is to convert the current structure of the CTD casts 
(from Andry Ratsimandresy) text files into yearly ones containing the following variables; time, 
latitude, longitude, temperature, salinity, level, instrument ID, and sounder depth

The input and output paths are local to the source computer, change as needed.

'''

##########################################################################################################
##(1) - Isolate and reformat casts according to year

#Declare the path to the files
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CTD_Andry/data_raw/'
path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CTD_Andry/data_processed/'

#Determine all the names of the folders
folders = np.array(os.listdir(path))
folders = folders[~(folders == 'individual_nc')]

#Declare the maximum depth for the .nc files
max_depth = 5000

#Declare the years covered
years = np.arange(2009,2020+1)

#Cycle through each folder
for folder in folders:

	#Determine the names of all files within the folder
	files = os.listdir(path + folder)

	#Cycle through each file
	for file in files:

		#Open the individual text file
		data = pd.read_csv(path + folder +'/'+ file, header=0, delimiter=r"\s+").to_dict('list')

		#Re-format the necessary variables 

		#Re-format the 1D variables
		time = np.unique(data['DeployTime(ISO8601)'])
		time = np.array([np.datetime64(i) for i in time])

		#Ensure that the events are used to isolate all 1D variables
		events = np.unique(data['Event'])
		time = [data['DeployTime(ISO8601)'][np.where(data['Event'] == i)[0][0]] for i in events]
		time_data = np.array([np.datetime64(i) for i in time])
		latitude = np.array([data['Latitude(DecDeg)'][np.where(data['Event'] == i)[0][0]] for i in events])
		longitude = np.array([data['Longitude(DecDeg)'][np.where(data['Event'] == i)[0][0]] for i in events])
		Cruise = np.array([data['Cruise'][np.where(data['Event'] == i)[0][0]] for i in events])
		Station = np.array([data['Station'][np.where(data['Event'] == i)[0][0]] for i in events])
		Type = np.array([data['Type'][np.where(data['Event'] == i)[0][0]] for i in events])
		BotDepth = np.array([data['BotDepth(m)'][np.where(data['Event'] == i)[0][0]] for i in events])
		Instrument = np.array([data['Instrument'][np.where(data['Event'] == i)[0][0]] for i in events])

		#Re-format the 2D variables
		temp_2D = np.zeros((events.size,max_depth))
		saln_2D = np.zeros((events.size,max_depth))

		#Cycle through each of the casts
		for i,value in enumerate(events):

			#Slice the temperature, pressure, and salinity
			#If pressure is present
			if np.isin('Pressure(dbar)',list(data)):
				pres_slice = np.array(data['Pressure(dbar)'])[np.array(data['Event']) == value]
			elif np.isin('Depth(m)',list(data)):
				pres_slice = np.array(data['Depth(m)'])[np.array(data['Event']) == value]

			temp_slice = np.array(data['Temperature(oC)'])[np.array(data['Event']) == value]
			saln_slice = np.array(data['Salinity'])[np.array(data['Event']) == value]

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

		#Remove the -99 values
		temp_2D[temp_2D <= -2] = np.nan
		saln_2D[saln_2D <= 0] = np.nan



		#Create individual .nc files for each of the text files
		#Set up the individual .nc file
		nc_out = nc.Dataset(path+'individual_nc/'+file[:-4]+'.nc','w')
	
		#File information
		nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
		nc_out.title = 'CTD Casts from Andry, Netcdf File' #Temporary title for the .nc file
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
		file_names = nc_out.createVariable('file_names', str, ('time'), zlib=True)
		Cruises = nc_out.createVariable('cruise', str, ('time'), zlib=True)
		Stations = nc_out.createVariable('station', str, ('time'), zlib=True)
		Types = nc_out.createVariable('type', str, ('time'), zlib=True)
		BotDepths = nc_out.createVariable('bottom_depth', np.float32, ('time'), zlib=True)
		Instruments = nc_out.createVariable('instrument', str, ('time'), zlib=True)

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
		latitudes[:] = latitude
		longitudes[:] = longitude
		times[:] = time_data
		file_names[:] = np.tile(file,time_data.size)
		Cruises[:] = Cruise.astype(str)
		Stations[:] = Station.astype(str)
		Types[:] = Type.astype(str)
		BotDepths[:] = BotDepth
		Instruments[:] = Instrument.astype(str)

		#Fill 2D structure
		temp[:,:] = temp_2D
		saln[:,:] = saln_2D

		#Convert to time stamps
		time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_data]
		times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
		levels[:] = np.arange(max_depth)

		#Save and close the .nc file
		nc_out.close()

		print(file+' done.')


#Now go about joining all like-year nc files together
files = os.listdir(path+'individual_nc/')

#Cycle through each of the files and open the ones with shared years
for year in years:

	nc_files = []
	#Cycle through each of the files
	for file in files:
		if file[:4] == str(year):
			nc_files.append(path+'individual_nc/'+file)

	#Merge all the files along the time axes
	nc_out = xr.open_mfdataset(nc_files,combine='nested',concat_dim = 'time')



	#Sort by time and save
	nc_out = nc_out.sortby('time')
	nc_out = nc_out.sel(time = nc_out['time.year'] == year)

	#Flatten the arrays to ensure they save properly
	for i in ['file_names','cruise','station','type','instrument']:
		nc_out[i] = nc_out[i].astype(str)

	nc_out.to_netcdf(path_output+str(year)+'.nc',encoding={'time':{'units': "seconds since 1900-01-01 00:00:00"}})
	nc_out.close()

	print(str(year)+' done.')
