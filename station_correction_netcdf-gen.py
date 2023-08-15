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








