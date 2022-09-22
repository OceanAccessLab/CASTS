import xarray as xr
import netCDF4 as nc
import numpy as np
from functools import reduce
import time as tt
import os
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

The purpose of this script is to combine the .nc filles, by year, from all sources
So far this includes only the BIO (pre and post 2010) and .pfiles sources
Script should be made in such a way that additional data can be easily added 

If a duplicate does exist in the netcdf_gen files, take the BIO_Climate data instead

'''
#######################################################################################


#Import the data by year
years = np.arange(1912,2021+1).astype(str)
path = {}
path['path1'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CCAS20/NetCDF/BIO_Climate_Databases/' #BIO, 1912-2010
path['path2'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CCAS20/NetCDF/database_2008-2017/' #BIO, 2008-2017
path['path3'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CCAS20/NetCDF/database_2018/' #BIO, 2018
path['path4'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/duplicates/' #NAFC, 1912-2021

#Set up the result destination
path_final = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined/'



#Rounding for the duplicate determination
time_rounder = 'datetime64[h]' #by hour
spatial_rounder = 2 #to 2 decimal points

#Cycle through each year
for year in years[:]:

	#Import each of the files, if one exists
	ds = {}

	#Cycle through each path
	for i in path.keys():
		if os.path.isfile(path[i]+str(year)+'.nc') == True:
			ds[i] = xr.open_dataset(path[i]+str(year)+'.nc')


	#If there is more than one dataset present for the year, proceed
	if np.size(list(ds.keys())) > 1:

		#Import the time, latitude, and longitude from each
		time,latitude,longitude = {},{},{}
		spatiotemporal_comp = {}
		duplicate_flag = {}

		#Cycle through each of the available paths
		for i in ds.keys():
			time[i] = ds[i].time.values.astype(time_rounder)
			latitude[i] = ds[i].latitude.values.round(spatial_rounder)
			longitude[i] = ds[i].longitude.values.round(spatial_rounder)

			#Create a spatio-temporal comp (in a string)
			spatiotemporal_comp[i] = np.array([
				time[i].astype(time_rounder).astype(str),
				latitude[i].round(spatial_rounder).astype(str),
				longitude[i].round(spatial_rounder).astype(str)
				])
			spatiotemporal_comp[i] = np.array([row[0]+','+row[1]+','+row[2] for row in spatiotemporal_comp[i].T])

			#Create empty flags array
			duplicate_flag[i] = np.full(time[i].size,False)

			#Determine where there are duplicates in each of the paths
			#There should be path1 for each year so don't compare like with like
			if (i == 'path1') == False:
				flags = np.isin(spatiotemporal_comp['path1'],spatiotemporal_comp[i])

				#Cycle through each of the duplicate locations 
				for ii in np.arange(flags.size):

					#Proceed if duplicate exists
					if flags[ii] == True:

						#Determine how many recordings are in path1
						place1 = (~np.isnan(ds['path1'].temperature[ii].values)).sum()

						#Determine how many recordings are in the orther path
						loc = np.where(spatiotemporal_comp[i] == spatiotemporal_comp['path1'][ii])[0][0]
						place2 = (~np.isnan(ds[i].temperature[loc].values)).sum()

						#If path1 has at least 90% of the recordings present in the other path, choose path1
						if (place1/ds[i].level.size+0.1) > place2/ds[i].level.size:
							
							#Mark the duplicates in path1 as False
							duplicate_flag['path1'][ii] = False
							duplicate_flag[i][loc] = True

						else:
							#Mark the duplicates in the other path as False
							duplicate_flag['path1'][ii] = True
							duplicate_flag[i][loc] = False
		
		#Compare the datasets that are not path1 with one another
		if np.size(list(ds.keys())) > 2:
			ds_keys_new = np.array([*ds])[1:]
			flags = np.isin(spatiotemporal_comp[ds_keys_new[0]],spatiotemporal_comp[ds_keys_new[1]])

			#Remove the flags where a True value is already present (removed)
			flags[duplicate_flag[ds_keys_new[0]]] = False

			#Cycle through each of the duplicate locations 
			for ii in np.arange(flags.size):

				#Proceed if duplicate exists
				if flags[ii] == True:

					#Determine how many recordings are in path1
					place1 = (~np.isnan(ds[ds_keys_new[0]].temperature[ii].values)).sum()

					#Determine how many recordings are in the orther path
					loc = np.where(spatiotemporal_comp[i] == spatiotemporal_comp[ds_keys_new[0]][ii])[0][0]
					place2 = (~np.isnan(ds[ds_keys_new[1]].temperature[loc].values)).sum()

					#Determine if either of the spots have been flagged already
					if duplicate_flag[ds_keys_new[0]][ii] == False and duplicate_flag[ds_keys_new[1]][loc] == False:

						#Take whichever path has more recordings
						if place1 > place2:
							
							#Mark the duplicates in path1 as False
							duplicate_flag[ds_keys_new[0]][ii] = False
							duplicate_flag[ds_keys_new[1]][loc] = True

						else:
							#Mark the duplicates in the other path as False
							duplicate_flag[ds_keys_new[0]][ii] = True
							duplicate_flag[ds_keys_new[1]][loc] = False


		#Remove the duplicates from each of the sources
		for i in ds.keys():
			ds[i] = ds[i].sel(time = ~duplicate_flag[i])

		#Some variables will need to be renamed for consistency
		for i in ds.keys():
			if i == 'path1':
				ds[i] = ds[i].rename({'datatype': 'instrument_ID',})
				ds[i] = ds[i].rename({'maximumdepth': 'sounder_depth'})		
			if i == 'path2':
				ds[i] = ds[i].rename({'stnid': 'instrument_ID'})
				ds[i] = ds[i].rename({'maximumdepth': 'sounder_depth'})						

		#Isolate for the variables of interest
		variables = [
		'level',
		'latitude',
		'longitude',
		'instrument_ID',
		'sounder_depth',
		'temperature',
		'salinity',
		]
		for i in ds.keys():
			ds[i] = ds[i][variables]

		#Determine which paths are still populated
		ds_keys_new = np.array([ds[i].time.size for i in ds.keys()]).astype(bool)
		ds_keys_new = np.array(list(ds.keys()))[ds_keys_new]

		#Merge the paths together into one dataset
		if ds_keys_new.size > 1:
			ds_merged = xr.concat([ds[i] for i in ds_keys_new],dim='time',combine_attrs='override')
		else:
			ds_merged = ds[ds_keys_new[0]]
		ds_merged = ds_merged.sortby('time')

		#Re-declare the file attributes
		ds_merged.attrs = {
		'Conventions': 'CF-1.6',
		'title': 'Canadian Atlantic Shelf Hydrography (CASH) Dataset',
		'institution': 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada',
		'source': 'https://github.com/OceanAccessLab/CASH',
		'references': 'Cyr et al., 2019',
		'description': 'Temperature and salinity cast records from the CASH Dataset.',
		'history': 'Created ' + tt.ctime(tt.time())
		}

		#Save the merged dataset
		ds_merged.to_netcdf(path_final+str(year)+'.nc',encoding={'time':{'units': "seconds since 1900-01-01 00:00:00"}})
		ds_merged.close()

	#If there is only one dataset present, copy then format and save
	else:

		#Rename variables depending on which path is present
		for i in ds.keys():
			if i == 'path1':
				ds[i] = ds[i].rename({'datatype': 'instrument_ID',})
				ds[i] = ds[i].rename({'maximumdepth': 'sounder_depth'})		
			if i == 'path2':
				ds[i] = ds[i].rename({'stnid': 'instrument_ID'})
				ds[i] = ds[i].rename({'maximumdepth': 'sounder_depth'})						

		#Isolate for the variables of interest
		variables = [
		'level',
		'latitude',
		'longitude',
		'instrument_ID',
		'sounder_depth',
		'temperature',
		'salinity',
		]
		for i in ds.keys():
			ds[i] = ds[i][variables]

		#Create the new ds
		ds_merged = ds[i]

		#Re-declare the file attributes
		ds_merged.attrs = {
		'Conventions': 'CF-1.6',
		'title': 'Canadian Atlantic Shelf Hydrography (CASH) Dataset',
		'institution': 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada',
		'source': 'https://github.com/OceanAccessLab/CASH',
		'references': 'Cyr et al., 2019',
		'description': 'Temperature and salinity cast records from the CASH Dataset.',
		'history': 'Created ' + tt.ctime(tt.time())
		}

		#Save the merged dataset
		ds_merged.to_netcdf(path_final+str(year)+'.nc',encoding={'time':{'units': "seconds since 1900-01-01 00:00:00"}})
		ds_merged.close()

	print(year+' done.')















time = ds_merged.time.values.astype(time_rounder)
latitude = ds_merged.latitude.values.round(spatial_rounder)
longitude = ds_merged.longitude.values.round(spatial_rounder)

#Create a spatio-temporal comp (in a string)
spatiotemporal_comp = np.array([
	time.astype(time_rounder).astype(str),
	latitude.round(spatial_rounder).astype(str),
	longitude.round(spatial_rounder).astype(str)
	])
spatiotemporal_comp = np.array([row[0]+','+row[1]+','+row[2] for row in spatiotemporal_comp.T])







