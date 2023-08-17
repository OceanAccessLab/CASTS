import xarray as xr
import netCDF4 as nc
import numpy as np
from functools import reduce
import time as tt
import os
import warnings
import pfile_tools as p
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt


#######################################################################################
'''	

The purpose of this script is to clean the station ID variable in the final netcdf 
files for each year.
Each station ID that will be included is in a .xlsx file. 
Rewrite this data if you want to include more stations.

'''
#######################################################################################


#Import the station ID data
station_path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Other/station_coordinates/'
station_ds = pd.read_excel(station_path+'Station_ID.xlsx')

#Isolate the data of interest
program = station_ds['Program'].values
station_ID = station_ds['Station_ID'].values
station_lat = station_ds['Latitude_decimal'].values
station_lon = station_ds['Longitude_decimal'].values*-1

#Define how much "padding" is given around station
padding = 0.0125

#Define which files will have the station IDs specified 
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined_region_filtered/'
years = np.arange(1912,2022+1).astype(str)

#Cycle through each year
for year in years:

	#Import the dataset
	ds = xr.open_dataset(path+year+'.nc')

	#Import the latitude, longitude, station ID
	ds_lat = ds.latitude.values 
	ds_lon = ds.longitude.values
	ds_ID = ds.station_ID.values

	#Define the new station_ID data
	new_ID = np.full(ds_ID.size,'',dtype='<U32')

	#Cycle through each station ID from the xlsx file
	for i in np.arange(station_ID.size):

		#Isolate which casts fall into them 
		place1 = (ds_lon >= (station_lon[i] - padding))*(ds_lon <= (station_lon[i] + padding))
		place2 = (ds_lat >= (station_lat[i] - padding))*(ds_lat <= (station_lat[i] + padding))

		#Write the new station ID
		new_ID[place1*place2] = program[i]+'_'+station_ID[i]

	#Record the new data
	ds['station_ID_manual'] = (('time'), new_ID)

	#LAST MINUTE CHANGE - THIS SHOULD BE MOVED TO FIRST STEP
	#Rename IEO 
	source = ds.source.values
	source[source == 'IEO-Spain'] = 'EU_NAFO'
	source[source == 'BIO-OMO'] = 'BIO-OMM'
	ds['source'][:] = source

	#Edit the attributes to be constant for each year
	#Time
	ds['time'].attrs['standard_name'] = 'time'
	ds['time'].attrs['variable_name'] = 'time'

	#Level
	ds['level'].attrs['standard_name'] = 'sea_water_pressure'
	ds['level'].attrs['maximum_level'] = '5000'
	ds['level'].attrs['units'] = 'dbar'
	ds['level'].attrs['variable_name'] = 'level'
	ds['level'].attrs['example_boundary'] = 'level 0: (0, 1], level 4900: (4900,5000]'

	#Temperature
	ds['temperature'].attrs['standard_name'] = 'sea_water_potential_temperature'
	ds['temperature'].attrs['variable_name'] = 'temperature'
	ds['temperature'].attrs['units'] = 'degree_C'
	ds['temperature'].attrs['valid_min'] = '-1.825'
	ds['temperature'].attrs['valid_max'] = '35'
	ds['temperature'].attrs['missing_value'] = 'nan'

	#Salinity
	ds['salinity'].attrs['standard_name'] = 'sea_water_practical_salinity'
	ds['salinity'].attrs['variable_name'] = 'salinity'
	ds['salinity'].attrs['units'] = 'psu'
	ds['salinity'].attrs['valid_min'] = '0'
	ds['salinity'].attrs['valid_max'] = '45'
	ds['salinity'].attrs['missing_value'] = 'nan'

	#Longitude
	ds['longitude'].attrs['standard_name'] = 'longitude'
	ds['longitude'].attrs['variable_name'] = 'longitude'
	ds['longitude'].attrs['units'] = 'degree_east'

	#Latitude
	ds['latitude'].attrs['standard_name'] = 'latitude'
	ds['latitude'].attrs['variable_name'] = 'latitude'
	ds['latitude'].attrs['units'] = 'degree_north'

	#Trip_ID
	ds['trip_ID'].attrs['standard_name'] = 'trip_identification'
	ds['trip_ID'].attrs['variable_name'] = 'trip_ID'

	#Source
	ds['source'].attrs['standard_name'] = 'cast_source'
	ds['source'].attrs['variable_name'] = 'source'

	#Instrument_ID
	ds['instrument_ID'].attrs['standard_name'] = 'instrument_identification'
	ds['instrument_ID'].attrs['variable_name'] = 'instrument_ID'

	#Instrument_type
	ds['instrument_type'].attrs['standard_name'] = 'instrument_type'
	ds['instrument_type'].attrs['variable_name'] = 'instrument_type'

	#File_names
	ds['file_names'].attrs['standard_name'] = 'file_names'
	ds['file_names'].attrs['variable_name'] = 'file_names'

	#Station_ID
	ds['station_ID'].attrs['standard_name'] = 'station_identification'
	ds['station_ID'].attrs['variable_name'] = 'station_ID'

	#Sounder_depth
	ds['sounder_depth'].attrs['standard_name'] = 'sounder_depth'
	ds['sounder_depth'].attrs['variable_name'] = 'sounder_depth'
	ds['sounder_depth'].attrs['units'] = 'dbar'

	#Station_ID_manual
	ds['station_ID_manual'].attrs['standard_name'] = 'station_identification_manual'
	ds['station_ID_manual'].attrs['variable_name'] = 'station_ID'

	#Global Attributes
	ds.attrs['title'] = 'Canadian Atlantic Shelf Temperature-Salinity (CASTS) Data Product, 2022, V1.0'
	ds.attrs['institution'] = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada (DFO)'
	ds.attrs['source'] = 'https://github.com/OceanAccessLab/CASH'
	ds.attrs['reference'] = 'Currently NA'
	ds.attrs['description'] = 'Temperature and salinity cast records from the CASTS Dataset'
	ds.attrs['comment'] = 'No constraints on data access or use'
	ds.attrs['Conventions'] = 'CF-1.6'
	ds.attrs['DOI'] = 'https://doi.org/10.20383/102.0739'
	ds.attrs['creator_names'] = 'Frederic Cyr, Jonathan Coyne'
	ds.attrs['creator_emails'] = 'frederic.cyr@dfo-mpo.gc.ca, jonathan.coyne@dfo-mpo.gc.ca'
	ds.attrs['geospatial_lon_max'] = '-42degE'
	ds.attrs['geospatial_lon_min'] = '-100degE'
	ds.attrs['geospatial_lat_min'] = '35degN'
	ds.attrs['geospatial_lat_max'] = '80degN'
	ds.attrs['file_year'] = year
	ds.attrs['keywords'] = 'Ocean Science, Ocean Temperature, Ocean Salinity, Historical Data'


	#Flatten the arrays to ensure they save properly
	for i in ['trip_ID','source','instrument_ID','instrument_type','file_names','station_ID','station_ID_manual','sounder_depth']:
		ds[i] = ds[i].astype(str)

	#Save the isolated netcdf files
	ds.to_netcdf('/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/final_product/'+\
		str(year)+'.nc',
		encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	ds.close()
	print(year+' done.')

#Cycle through all the files and compress them
#Right now this needs to be done from nc_old after files have been moved there manually, not cool
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/final_product/'
path_final = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/final_product/compressed/'
for year in np.arange(1912,2022+1).astype(str)[:]:

	#Compress all files for memory purposes
	exp = 'nccopy -d 9 '+path+year+'.nc '+path_final+year+'.nc'
	os.system(exp)

	print(year+' done.')








