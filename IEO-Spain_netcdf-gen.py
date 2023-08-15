import xarray as xr
import netCDF4 as nc
import numpy as np
from functools import reduce
import time as tt
import os
import csv
import pfile_tools as p
import pandas as pd
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt


'''

The purpose of this script is convert the IEO-Spain CTD data from the original 
.int file type to a uniform yearly formatted NetCDF file type.

The input and output paths are specific to the source computer, change as needed

'''

##########################################################################################################
##(1) - Convert each of the folder files to a netcdf

#Define the path of data
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/IEO_Spain/data_raw/'
path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/IEO_Spain/data_processed/'
max_depth = 5000

#Determine all the file names (end with .dat)
files = np.array(os.listdir(path))
filt = np.zeros(files.size).astype(bool)
for i,value in enumerate(files):
	if value[-3:] == 'dat':
		filt[i] = True
files = files[filt]

#Cycle through each of the files
for file in files:

	#Load in the text file
	with open(path+file, encoding='ISO-8859-1') as line:
		reader = csv.reader(line, delimiter=':')
		cast_head = [row for row in reader]
	cast_head = np.array(cast_head,dtype=object)

	#The text '*SURFACE SAMPLES=' defines the beginning of a cast
	#The recording -999 for PRES marks the end of the cast
	#In between is metadata
	cast_start = []
	cast_end = []
	for i,value in enumerate(cast_head):
		if value == ['*SURFACE SAMPLES=']:
			cast_start.append(i+2)

			#Determine the end of the cast
			for ii in np.arange(i+3,cast_head.size):

				if np.size(cast_head[ii][0].split()) != np.size(cast_head[i+3][0].split()):
					cast_end.append(ii)
					break
	cast_end.append(cast_head.size)

	#Determine the location of all metadata 
	header_start = []
	header_end = []
	header_start.append(0)
	for i in cast_end[:-1]:
		header_start.append(i+1)
	for i in cast_start:
		header_end.append(i-2)

	#Cycle through each of the identified casts
	#Set up all the 1D variables
	time_1D = np.zeros(np.size(cast_start)).astype('datetime64[m]')
	lat_1D = np.zeros(np.size(cast_start)).astype(float)
	lon_1D = np.zeros(np.size(cast_start)).astype(float)
	depth_1D = np.zeros(np.size(cast_start)).astype(str)
	QC_1D = np.zeros(np.size(cast_start)).astype(str)
	instrument_id_1D = np.zeros(np.size(cast_start)).astype(str)
	datatype_1D = np.zeros(np.size(cast_start)).astype(str)
	station_id_1D = np.zeros(np.size(cast_start)).astype(str)

	#Set up all the 2D variables
	temp_2D = np.full((np.size(cast_start),max_depth), np.nan)
	saln_2D = np.full((np.size(cast_start),max_depth), np.nan)


	#Populate the 1D variables 
	for i in np.arange(np.size(header_start)):

		#Cycle through each of the header lines
		for ii in np.arange(header_start[i],header_end[i]):

			#Isolate time, latitude, longitude, depth, QC
			if cast_head[ii][0].startswith('*DATE='):
				place1 = cast_head[ii][0].split()

				#Time
				minute = place1[1][-2:]
				hour = place1[1][-4:-2]
				day = place1[0][-8:-6]
				month = place1[0][-6:-4]
				year = place1[0][-4:]
				time_1D[i] = np.datetime64(year+'-'+month+'-'+day+'T'+hour+':'+minute)

				#Latitude and Longitude
				lat_1D[i] = float(place1[2][5:]) + (float(place1[3])/60.0)
				if place1[2][4] == 'S':
					lat_1D[i] = lat_1D[i]*-1
				lon_1D[i] = float(place1[4][5:]) + (float(place1[5])/60.0)
				if place1[4][4] == 'W':
					lon_1D[i] = lon_1D[i]*-1

				#Depth 
				depth_1D[i] = place1[6][6:]

				#Quality Control
				QC_1D[i] = place1[7][3:]

			#Isolate instrument id, data type
			if cast_head[ii][0].startswith('*DC HISTORY='):
				place1 = cast_head[ii][0][12:].split()

				#Instrument id
				if np.size(place1) >= 4:
					instrument_id_1D[i] = place1[0]+'-'+place1[1]
				else:
					instrument_id_1D[i] = ''

				#Data Type
				if np.size(place1) >= 4:
					datatype_1D[i] = place1[2]
				else:
					#We know this data is CTD
					datatype_1D[i] = 'CTD'

			#Isolate the station id 
			if cast_head[ii][0].startswith('*DM HISTORY='):
				if np.size(cast_head[ii]) >= 2:
					station_id_1D[i] = cast_head[ii][1]
				else:
					station_id_1D[i] = ''

	#Populate the 2D variables
	for i in np.arange(np.size(header_start)):

		#Cycle through each of the casts
		for ii in np.arange(cast_start[i]+1,cast_end[i]-1): 
			place1 = cast_head[ii][0].split()

			#Fill in pressure, temperature, and salinity
			temp_2D[i,int(place1[0])] = float(place1[1])
			saln_2D[i,int(place1[0])] = float(place1[2])


	#Save in NetCDF format
	#Set up the .nc file
	nc_out = nc.Dataset(path_output+file[:-4]+'.nc','w')
	
	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'IEO-Spain Yearly Netcdf File' #Temporary title for the .nc file
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
	depths = nc_out.createVariable('depth', str, ('time'), zlib=True)
	QCs = nc_out.createVariable('QC', str, ('time'), zlib=True)
	instrument_ids = nc_out.createVariable('instrument_id', str, ('time'), zlib=True)
	datatypes = nc_out.createVariable('datatype', str, ('time'), zlib=True)
	station_ids = nc_out.createVariable('station_id', str, ('time'), zlib=True)
	file_names = nc_out.createVariable('file_names', str, ('time'), zlib=True)
	
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

	#Sort according to time
	filt = np.argsort(time_1D)

	#Fill in the 1D variables 
	latitudes[:] = lat_1D[filt]
	longitudes[:] = lon_1D[filt]
	times[:] = time_1D[filt]
	depths[:] = depth_1D[filt]
	QCs[:] = QC_1D[filt]
	instrument_ids[:] = instrument_id_1D[filt]
	datatypes[:] = datatype_1D[filt]
	station_ids[:] = station_id_1D[filt]
	file_names[:] = np.tile(file, time_1D.size)

	#Fill 2D structure
	temp[:,:] = temp_2D[filt]
	saln[:,:] = saln_2D[filt]

	#Convert to time stamps
	time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_1D[filt]]
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_depth)

	#Save and close the .nc file
	nc_out.close()
	print(file+' done.')


##########################################################################################################
##(2) - Cycle through and isolate by year

#Now cycle through all the NetCDF files and save by year
ds = {}
for file in files:
	ds[file[:-4]] = xr.open_dataset(path_output+file[:-4]+'.nc')

#Merge all of the sources together
ds_merged = xr.concat([ds[file[:-4]] for file in files],dim='time',combine_attrs='override')
ds_merged = ds_merged.sortby('time')

#Determine which years are covered
years = np.unique(ds_merged['time.year'])

#Cycle through each year
for year in years:

	#Slice the ds_merged by year
	ds_isolate = ds_merged.isel(time=ds_merged['time.year'] == year)
	ds_isolate.attrs = {
	'Conventions': 'CF-1.6',
	'title': 'IEO-Spain Yearly NetCDF File',
	'institution': 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada',
	'source': 'https://github.com/OceanAccessLab/CASH',
	'references': 'Cyr et al., 2022',
	'description': 'Temperature and salinity cast records from the CASTS Dataset.',
	'history': 'Created ' + tt.ctime(tt.time())
	}

	#Save the merged dataset
	path_final = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/IEO_Spain/netcdf_yearly/'
	ds_isolate.to_netcdf(path_final+str(year)+'.nc',\
		encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	ds_isolate.close()
	print(str(year)+' done.')





