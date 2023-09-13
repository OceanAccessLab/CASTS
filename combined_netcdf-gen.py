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

The purpose of this script is to combine the .nc filles, by year, from all sources
So far this includes only the BIO (pre and post 2010) and .pfiles sources
Script should be made in such a way that additional data can be easily added 

If a duplicate does exist in the netcdf_gen files, take the BIO_Climate data instead

'''
#######################################################################################


#Import the data by year
years = np.arange(1912,2023+1).astype(str)
path = {}
path['path1'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CCAS20/NetCDF/BIO_Climate_Databases/' #BIO, 1913-2010
path['path2'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CCAS20/NetCDF/database_2008-2017/' #BIO, 2008-2017
path['path3'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CCAS20/NetCDF/database_2018/' #BIO, 2018
path['path4'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/duplicates/' #NAFC, 1999-2022
path['path5'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CIOOS_ERRDAP/data_processed/' #CIOOS-ERRDAP, 1996-2020
path['path6'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NEFSC/yearly_netcdf_files/' #NEFSC, 1981-2021
path['path7'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CTD_Andry/data_processed/' #CTD Andry, 2009-2020
path['path8'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/ASCII/data_processed/' # MLI-Sourced, 1972-2022
path['path9'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NCEI_GTSPP/data_raw/data_processed/bottom_flagged/' # NCEI-GTSPP, 1990-2019
path['path10'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Polar_Data_Catalogue/data_product/netcdf_yearly/' # Polar Data Catalogue, 2002-2020
path['path11'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/IEO_Spain/netcdf_yearly/' # IEO Spain, 2019-2021
path['path12'] = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Other/Freds_Files/netcdf_yearly/' #Marine Institute, 2000-2008

#Declare the path sources
path_source = {}
path_source['path1'] = 'Climate'
path_source['path2'] = 'BIO-OMM'
path_source['path3'] = 'BIO-OMM'
path_source['path4'] = 'NAFC-Oceanography'
path_source['path5'] = 'CIOOS'
path_source['path6'] = 'NEFSC'
path_source['path7'] = 'NAFC-Aquaculture'
path_source['path8'] = 'MLI'
path_source['path9'] = 'NCEI'
path_source['path10'] = 'Polar-Data-Catalogue'
path_source['path11'] = 'EU-NAFO'
path_source['path12'] = 'Marine-Institute-NL'

#Set up the result destination
path_final = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined/'

#Rounding for the duplicate determination
time_rounder = 'datetime64[m]' #by minute
spatial_rounder = 2 #to 2 decimal points



#Cycle through each year
for year in years[:]:

	#Import each of the files, if one exists
	ds = {}

	#Cycle through each path
	for i in path.keys():
		if os.path.isfile(path[i]+str(year)+'.nc') == True:
			ds[i] = xr.open_dataset(path[i]+str(year)+'.nc')

	#Some variables will need to be renamed for consistency
	for i in ds.keys():
		if i == 'path1':
			ds[i] = ds[i].rename({'datatype': 'instrument_ID',})
			ds[i] = ds[i].rename({'maximumdepth': 'sounder_depth'})		
		if i == 'path2' or i == 'path3':
			ds[i] = ds[i].rename({'maximumdepth': 'sounder_depth'})		
		if i == 'path3':
			ds[i] = ds[i].rename({'datatype': 'instrument_ID'})
		if i == 'path6':
			ds[i] = ds[i].rename({'CTD_serial_number': 'instrument_ID'})
			ds[i] = ds[i].rename({'deployment_method': 'instrument_type'})
			ds[i] = ds[i].rename({'station_number': 'comments'})
			ds[i]['trip_ID'] = (['time'], 
					np.char.array(ds[i].country_code.values.astype(str))+'_'+\
					np.char.array(ds[i].ship_code.values.astype(str))+'_'+\
					np.char.array(ds[i].cruise_code.values.astype(str))
					)
		if i == 'path9':
			ds[i] = ds[i].rename({'file_name': 'file_names'})
			ds[i] = ds[i].rename({'source_ID': 'comments'})
			ds[i] = ds[i].rename({'data_type': 'instrument_ID',})
		if i == 'path11':
			ds[i] = ds[i].rename({'datatype': 'instrument_ID'})
			ds[i] = ds[i].rename({'instrument_id': 'instrument_type'})
			ds[i] = ds[i].rename({'station_id': 'comments'})
		if i == 'path12':
			ds[i] = ds[i].rename({'mission_id': 'trip_ID'})


	#Isolate for the variables of interest
	variables = np.array([
		'level',
		'latitude',
		'longitude',
		'instrument_ID',
		'instrument_type',
		'comments',
		'trip_ID',
		'file_names',
		'sounder_depth',
		'temperature',
		'salinity',
		])

	for i in ds.keys():

		#Check to see if the variable is present
		variables_temp = variables[np.isin(variables,list(ds[i].variables))]
		ds[i] = ds[i][variables_temp]

		#If the variable is not available make a nan one
		if np.isin(variables,list(ds[i].variables)).all() == False:
			for ii in variables[~np.isin(variables,list(ds[i].variables))]:
				ds[i][ii] = (['time'], np.full(ds[i].time.size, np.nan).astype(str))


	#Merge all of the variables together
	ds_merged = xr.concat([ds[i] for i in ds.keys()],dim='time',combine_attrs='override')
	ds_merged = ds_merged.rename({'comments': 'station_ID'})

	#Determine where each path is located
	path_location = np.concatenate([np.full(ds[i].time.size,path_source[i]) for i in ds.keys()])
	ds_merged['source'] = (['time'],path_location.astype(object))

	#Remove all casts which have a standard deviation of 0 throughout, if they have more than 2 measurements
	#Determine for both temperature and salinity
	warnings.simplefilter("ignore", category=RuntimeWarning)
	temp_std = ds_merged.temperature.std(axis=1).values == 0
	temp_nom = ds_merged.level.size - np.isnan(ds_merged['temperature']).sum(axis=1)
	saln_std = ds_merged.salinity.std(axis=1).values == 0
	saln_nom = ds_merged.level.size - np.isnan(ds_merged['salinity']).sum(axis=1)

	#Remove the constants
	ds_merged['temperature'][temp_std*(temp_nom > 3),:] = np.nan
	ds_merged['salinity'][saln_std*(saln_nom > 3),:] = np.nan

	#Remove the new empties
	temp_empties = np.isnan(ds_merged['temperature']).sum(axis=1) == ds_merged.level.size
	saln_empties = np.isnan(ds_merged['salinity']).sum(axis=1) == ds_merged.level.size
	ds_merged = ds_merged.sel(time = ~(temp_empties.values*saln_empties.values))

	#Define the time and distance radius around the point 
	time_cutoff = np.timedelta64(10,'m')
	dist_cutoff = 0.0125

	#Import the time, latitude, and longitude
	time = ds_merged.time.values
	latitude = ds_merged.latitude.values
	longitude = ds_merged.longitude.values

	#Create a variable marking duplicates
	duplicates_flag = np.full(ds_merged.time.shape,False) 

	#Cycle through each point 
	for i in np.arange(time.size):

		#Determine if the cast is already flagged as a duplicate, if so, skip
		if duplicates_flag[i] == False:

			#Determine the distance and time between that point and every other
			time_diff = time - time[i]
			lon_diff = longitude - longitude[i]
			lat_diff = latitude - latitude[i]

			#Determine which points make the cutoffs
			time_diff = np.abs(time_diff) < time_cutoff
			lon_diff = np.abs(lon_diff) < dist_cutoff
			lat_diff = np.abs(lat_diff) < dist_cutoff

			#Determine where all three criteras are met
			flagged_casts = time_diff*lat_diff*lon_diff

			#Determine if there are duplicates present, if so proceed
			if flagged_casts.sum() > 1:

				#Determine the corresponding sources, index
				flagged_source = ds_merged.source[flagged_casts].values
				flagged_index = np.where(flagged_casts == True)[0]

				#Determine the number of temperature and salinity measurements in each
				temp = ds_merged.temperature[flagged_casts].values
				saln = ds_merged.salinity[flagged_casts].values

				#Determine the number of recordings in each cast
				nom_temp = (~np.isnan(temp)).sum(axis=1)
				nom_saln = (~np.isnan(saln)).sum(axis=1)

				#If Climate is present
				if np.isin('Climate',flagged_source) == True:

					#If Climate is present once
					if np.isin(flagged_source,'Climate').sum() == 1:
					
						#If the Climate location has at least 90% of recordings, keep
						clim_place1 = {}
						place1 = {}

						#Determine the number of recordings in all  spots
						for ii,value in enumerate(flagged_source):
							if value != 'Climate':
								place1[ii] = int(nom_temp[ii] + nom_saln[ii])
							else:
								clim_place1[value] = int(nom_temp[ii] + nom_saln[ii])

						#Determine if Climate has enough to be kept (90% of others)
						place1_check = np.zeros(flagged_casts.sum()-1).astype(bool)
						x = 0
						for ii,value in enumerate(flagged_source):
							if value != 'Climate':

								#Check number of measurements
								if place1[ii]*0.9 > clim_place1['Climate']:
									place1_check[x] = True
								x += 1

						#If all aren't, keep Climate, discard others
						if not any(place1_check):
							duplicates_flag[flagged_index[~(flagged_source == 'Climate')]] = True

						#Else take the cast with the highest number
						else:
							duplicates_flag[flagged_index[flagged_index != flagged_index[np.argmax(nom_temp+nom_saln)]]] = True

					#If Climate is the only source present
					elif np.isin(flagged_source,'Climate').all():
						duplicates_flag[flagged_index[flagged_index != flagged_index[np.argmax(nom_temp+nom_saln)]]] = True

					#If Climate is present more than twice, but other variables are present
					else:

						#If the Climate location has at least 90% of recordings, keep
						clim_place1 = {}
						place1 = {}
						
						#Determine the number of recordings in all  spots
						for ii,value in enumerate(flagged_source):
							if value != 'Climate':
								place1[ii] = int(nom_temp[ii] + nom_saln[ii])
							else:
								clim_place1[ii] = int(nom_temp[ii] + nom_saln[ii])

						#Determine which group has the highest number
						place1_winner = np.max([place1[ii]*0.9 for ii in place1])
						clim_place1_winner = np.max([clim_place1[ii] for ii in clim_place1])

						#If place1 group is higher
						if place1_winner > clim_place1_winner:
							duplicates_flag[flagged_index[flagged_index != flagged_index[np.argmax(nom_temp+nom_saln)]]] = True

						#If clim_place1 group is higher
						else:
							place1 = np.argmax((nom_temp+nom_saln)[flagged_source == 'Climate'])
							#Remove the non-Climate casts
							duplicates_flag[flagged_index[np.isin(flagged_index,flagged_index[flagged_source == 'Climate'],invert=True)]] = True
							#Remove the Climate casts that didn't win
							duplicates_flag[flagged_index[flagged_source == 'Climate'][flagged_index[flagged_source == 'Climate'] != flagged_index[flagged_source == 'Climate'][place1]]] = True

				#If Climate is not present
				else:
					#Determine which cast has the highest number of measurements and keep that one
					duplicates_flag[flagged_index[flagged_index != flagged_index[np.argmax(nom_temp+nom_saln)]]] = True


	#Remove the duplicates from the ds
	ds_merged = ds_merged.sel(time=~duplicates_flag)
	ds_merged = ds_merged.sortby('time')

	#Manually set stations of interest
	stns_of_interest = {
	'S27-01':	[47.54667,-52.58667],
	'HL-02':	[44.2670,-63.3170],
	'Prince-5':	[44.9300,-66.8500],
	}
	#Set the size box around which to isolate 
	padding = 0.0125

	#Cycle through each station
	for i in stns_of_interest:

		#Define the box coordinates
		lon_min = stns_of_interest[i][1]-padding 
		lon_max = stns_of_interest[i][1]+padding
		lat_min = stns_of_interest[i][0]-padding
		lat_max = stns_of_interest[i][0]+padding

		#Isolate the latitudes and longitudes of interest
		lat_flag = (ds_merged.latitude.values > lat_min)*(ds_merged.latitude.values < lat_max)
		lon_flag = (ds_merged.longitude.values > lon_min)*(ds_merged.longitude.values < lon_max)

		#Change the comments value for flagged casts
		ds_merged['station_ID'][lat_flag*lon_flag] = i

	#Check to ensure all points are within the desired year
	ds_merged = ds_merged.sel(time = ds_merged['time.year'] == int(year))

	#Re-declare the file attributes
	ds_merged.attrs = {
	'Conventions': 'CF-1.6',
	'title': 'Canadian Atlantic Shelf Temperature-Salinity (CASTS) Dataset',
	'institution': 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada',
	'source': 'https://github.com/OceanAccessLab/CASH',
	'references': 'Cyr et al., 2022',
	'description': 'Temperature and salinity cast records from the CASTS Dataset.',
	'history': 'Created ' + tt.ctime(tt.time())
	}

	#Flatten the arrays to ensure they save properly
	'''
	for i in ['trip_ID','source','instrument_ID','instrument_type','file_names','station_ID','sounder_depth']:
		ds_merged[i] = ds_merged[i].astype(str)
	'''

	#Save the merged dataset
	ds_merged.to_netcdf(path_final+str(year)+'.nc',\
		encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	ds_merged.close()
	print(year+' done.')


#Next, go through and isolate for the region of interest
#35 to 80 degrees North and 42 to 100 degrees West

#Cycle through each year
for year in years:

	#Open the netcdf file
	ds = xr.open_dataset(path_final+str(year)+'.nc')

	#Perform the cutoff for the temperature and salinity
	#Pre-define cut-offs for each variable that are non-realistic
	cutoff = {
	'temperature':	[-1.825,35],
	'salinity':		[0,45],
	}

	#Isolate the data
	place1 = ds.where(ds['temperature'] >= cutoff['temperature'][0], drop=False)  
	place1 = place1.where(place1['temperature'] <= cutoff['temperature'][1], drop=False)  
	place2 = ds.where(ds['salinity'] >= cutoff['salinity'][0], drop=False)  
	place2 = place2.where(place2['salinity'] <= cutoff['salinity'][1], drop=False)  
	ds['temperature'] = place1.temperature
	ds['salinity'] = place2.salinity

	#Isolate for the region of interest
	ds = ds.sel(time = (ds.latitude >= 35)*(ds.latitude <= 80))
	ds = ds.sel(time = (ds.longitude >= -100)*(ds.longitude <= -42))

	#TEMPORARY
	#Remove NAFC-Oceanography casts startiong with 2017_79 during 2017 or 2012_39 for 2012
	if np.isin(year, ['2017','2018']):

		#Write down the file name
		file_names = ds.file_names.values
		place1 = np.zeros(file_names.size).astype(bool)
		for i,value in enumerate(file_names):
			if str(value).startswith('2017_79'):
				place1[i] = True

		#Remove the flagged casts
		ds = ds.sel(time = (~place1))


	#Bin-average the values over specified depths
	bin1 = np.arange(0,1000+1,1)
	bin2 = np.arange(1010,2000+10,10)
	bin3 = np.arange(2100,5000+100,100)
	bins = np.concatenate((bin1,bin2,bin3))

	#Bin-average the depths
	temp_bins = ds.temperature.groupby_bins('level',bins=bins,include_lowest=True).mean()
	saln_bins = ds.salinity.groupby_bins('level',bins=bins,include_lowest=True).mean()

	ds_new = temp_bins.to_dataset()
	ds_new['salinity'] = saln_bins
	for i in [
		'longitude',
		'latitude',
		'trip_ID',
		'source',
		'instrument_ID',
		'instrument_type',
		'file_names',
		'station_ID',
		'sounder_depth']:
		ds_new[i] = ds[i]
	ds_new.attrs = ds.attrs
	ds_new = ds_new.assign_coords({'level_bins':('level_bins',bins[:-1])})
	ds_new = ds_new.rename({'level_bins':'level'})


	#Flatten the arrays to ensure they save properly
	for i in ['trip_ID','source','instrument_ID','instrument_type','file_names','station_ID','sounder_depth']:
		ds_new[i] = ds_new[i].astype(str)

	#Save the isolated netcdf files
	ds_new.to_netcdf('/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined_region/'+\
		str(year)+'.nc',
		encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	ds_new.close()
	print(year+' done.')
































