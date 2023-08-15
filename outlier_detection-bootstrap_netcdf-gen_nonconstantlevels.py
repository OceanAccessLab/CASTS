import xarray as xr
import netCDF4 as nc
import numpy as np
import time as tt
import warnings
import os
import math
import random
import pfile_tools as p
import pandas as pd
from scipy import stats,interpolate
from scipy.signal import butter, lfilter, freqz
import scipy
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt


#######################################################################################
'''

Determine data records which are outliers (both spatially and temporally)
The spatial factor needs to consider the distance between data points

There will be some human-decided factors that need to be kept track of for this method

This version of the code is built to handle the non-constant vertical spacing of 
newer CASTS versions.

'''
#######################################################################################


#######################################################################################
#(1)	-	Create Climatology

#Meshgrid lat/lon
#Define lat and lon on a constant grid across years
#According to the BIO_CLIMATE data paper latmin=35,latmax=80,lonmin=-100,lonmax=-42
dc = 1
x = np.arange(-100, -42, dc)
y = np.arange(35, 80, dc)  
lon, lat = np.meshgrid(x,y)

#Pre-define through each year
years = np.arange(1912,2022+1)

#Import the netcdf files
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined_depth_filtered/'
ds = xr.open_mfdataset(path+'*.nc')

#Load in the latitude, longitude, depth, time, temperature, salinity
ds_time = ds['time.month'].values
ds_lons = ds.longitude.values
ds_lats = ds.latitude.values
ds_temp = ds.temperature.values
ds_saln = ds.salinity.values


#Cycle through each month
for i in np.arange(12)[:]:

	#Create a monthly filter
	month_filt = ds_time == i+1

	#Create an empty array for all the binned measurements
	level = ds.level.values
	level_size = ds.level.size
	temp_3D_mean = np.full((level_size,y.shape[0]-1,x.shape[0]-1),np.nan)
	saln_3D_mean = np.full((level_size,y.shape[0]-1,x.shape[0]-1),np.nan)
	temp_3D_stdv = np.full((level_size,y.shape[0]-1,x.shape[0]-1),np.nan)
	saln_3D_stdv = np.full((level_size,y.shape[0]-1,x.shape[0]-1),np.nan)

	#Determine where each cast occurs throughout the month
	lats_slice = ds_lats[month_filt]
	lons_slice = ds_lons[month_filt]
	temp_month = ds_temp[month_filt]
	saln_month = ds_saln[month_filt]

	#Determine where each cast is 
	place2 = stats.binned_statistic_2d(
		lons_slice,
		lats_slice,
		temp_month[:,0],
		bins=[lon[0,:],lat[:,0]],
		statistic='count',expand_binnumbers=True
		)

	#Cycle through each point
	for ii in np.arange(place2.statistic.T.shape[0]):
		for iii in np.arange(place2.statistic.T.shape[1]):

			#Determine if there are more than three measurements
			if place2.statistic.T[ii,iii] > 3:

				#Isolate for the relevant casts
				cast_isolate = (place2.binnumber[0] == iii+1)*(place2.binnumber[1] == ii+1)	

				#Define the number of cycles
				noc = 100 #number of cycles
				nom = 100 #number of measurements

				#Define the depth window
				depth_window = 10

				#Cycle through each depth
				for iv in np.arange(level_size):

					#First, determine the vertical spacing
					if np.where(level >= level[iv] + depth_window)[0].size >= 1:
						below_spacing = np.where(level >= level[iv] + depth_window)[0][0]
						if level[below_spacing] - level[iv] > depth_window:
							below_spacing = iv
					else:
						below_spacing = iv
					if np.where(level <= level[iv] - depth_window)[0].size >= 1:
						above_spacing = np.where(level <= level[iv] - depth_window)[0][-1]
						if level[iv] - level[above_spacing] > depth_window:
							above_spacing = iv
					else:
						above_spacing = 0

					#Next, determine the temp depth window
					temp_slice = temp_month[cast_isolate][:,above_spacing:below_spacing+1]
					warnings.simplefilter("ignore", category=RuntimeWarning)
					temp_slice = np.nanmean(temp_slice, axis=1)

					if temp_slice[~np.isnan(temp_slice)].size > 3:

						#Define an array where the means are stored
						temp_means = np.zeros(noc)

						#Cycle through each interation
						for v in np.arange(noc):

							#Create a record of the chosen temperature
							temp_means[v] = np.mean(random.choices(
								temp_slice[~np.isnan(temp_slice)],
								k=nom))

						#Record the mean and the standard deviation
						temp_3D_mean[iv,ii,iii] = temp_means.mean()
						temp_3D_stdv[iv,ii,iii] = temp_means.std()

					#Determine if there's enough measurements to proceed, salinity
					saln_slice = saln_month[cast_isolate][:,above_spacing:below_spacing+1]
					warnings.simplefilter("ignore", category=RuntimeWarning)
					saln_slice = np.nanmean(saln_slice, axis=1)


					if saln_slice[~np.isnan(saln_slice)].size > 3:

						#Define an array where the means are stored
						saln_means = np.zeros(noc)

						#Cycle through each interation
						for v in np.arange(noc):

							#Create a record of the chosen temperature
							saln_means[v] = np.random.choice(
								saln_slice[~np.isnan(saln_slice)],
								size=nom).mean()

						#Record the mean and the standard deviation
						saln_3D_mean[iv,ii,iii] = saln_means.mean()
						saln_3D_stdv[iv,ii,iii] = saln_means.std()
					
					print('Depth: '+str(level[iv])+', Y: '+str(ii)+', X: '+str(iii)+', Month: '+str(i+1))

	#Save the temperature and salinity for each month 
	#Save the monthly climatology mean and stdv in a .nc format
	path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/climatology/'

	#Write a new .nc file
	nc_out = nc.Dataset(path+'bootstrap_climatology_'+"%.2d" % (i+1)+'.nc','w')

	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'Gridded Mean/Stdv temperature and Salinity using bootstrap method, combined_region'  #Temporary title for the .nc file
	nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
	nc_out.source = 'https://github.com/jn533213/AZMP_Fred/blob/main/Trial2_RData.py'
	nc_out.references = 'Fill later'
	nc_out.description = 'A binned-temperature,salinity test by jonathan.coyne@dfo-mpo.gc.ca'
	nc_out.history = 'Created ' + tt.ctime(tt.time())

	#Create dimensions
	level = nc_out.createDimension('level', level_size) 
	longitude = nc_out.createDimension('longitude', x.size-1) 
	latitude = nc_out.createDimension('latitude', y.size-1) 

	#Create coordinate variables
	levels = nc_out.createVariable('level', np.int32, ('level',))
	latitudes = nc_out.createVariable('latitude', np.float64, ('latitude',))
	longitudes = nc_out.createVariable('longitude', np.float64, ('longitude',))

	#Create the temperature and salinity variables
	temp_mean = nc_out.createVariable(
		'temperature_mean', 
		np.float64, 
		('level','latitude','longitude'), 
		zlib=True, fill_value=-9999)
	temp_ster = nc_out.createVariable(
		'temperature_ster', 
		np.float64,   
		('level','latitude','longitude'), 
		zlib=True, fill_value=-9999)
	saln_mean = nc_out.createVariable(
		'salinity_mean', 
		np.float64, 
		('level','latitude','longitude'), 
		zlib=True, fill_value=-9999)
	saln_ster = nc_out.createVariable(
		'salinity_ster', 
		np.float64, 
		('level','latitude','longitude'), 
		zlib=True, fill_value=-9999)

	# Variable Attributes
	latitudes.units = 'degree_north'
	longitudes.units = 'degree_east'
	levels.units = 'dbar'
	levels.standard_name = "pressure"
	levels.valid_min = 0
	temp_mean.units = 'Celsius'
	temp_mean.long_name = "Water Temperature Bootstrap Mean"
	temp_mean.standard_name = "sea_water_temperature_mean"
	temp_ster.units = 'Celsius'
	temp_ster.long_name = 'Water Temperature Bootstrap Standard Error'
	temp_ster.standard_name = 'sea_water_temperature_standard_error'
	saln_mean.units = 'psu'
	saln_mean.long_name = 'Salinity Bootstrap Mean'
	saln_mean.standard_name = 'salinity_mean'
	saln_ster.units = 'psu'
	saln_ster.long_name = 'Salinity Bootstrap Standard Error'
	saln_ster.standard_name = 'salinity_standard_error'

	#Save the temperature and salinity
	temp_mean[:,:,:] = temp_3D_mean
	saln_mean[:,:,:] = saln_3D_mean
	temp_ster[:,:,:] = temp_3D_stdv
	saln_ster[:,:,:] = saln_3D_stdv

	#Convert to time stamps
	levels[:] = ds.level.values
	latitudes[:] = y[:-1]
	longitudes[:] = x[:-1]

	nc_out.close()




#######################################################################################
#(2)	-	Determine Outliers


#Cycle through the years
for year in np.arange(1912,2022+1).astype(str)[:]:

	#Import the netcdf dataset, one year at a time
	path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined_depth_filtered/'
	#path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/duplicates/'
	ds = xr.open_dataset(path+str(year)+'.nc')

	#Determine the coresponding coordinates for each cast in the climatology
	#Define the lat/lon for the NAFC
	NAFC_lat = ds.latitude.values
	NAFC_lon = ds.longitude.values

	#Import the lat/lon from the climatology
	ds_clim = xr.open_dataset('/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/climatology/bootstrap_climatology_01.nc')
	#ds_clim = xr.open_dataset('/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/climatology/bootstrap_climatology_01.nc')	
	clim_lon,clim_lat = np.meshgrid(ds_clim.longitude.values, ds_clim.latitude.values)
	indx_lon,indx_lat = np.meshgrid(np.arange(clim_lon.shape[1]),np.arange(clim_lat.shape[0]))

	#Determine the indices for each of the NAFC casts
	NAFC_xcoord = interpolate.griddata(
	np.array((clim_lon.flatten(),clim_lat.flatten())).T,
	indx_lon.flatten().T,
	(NAFC_lon,NAFC_lat),
	method='nearest'
	)
	NAFC_ycoord = interpolate.griddata(
	np.array((clim_lon.flatten(),clim_lat.flatten())).T,
	indx_lat.flatten().T,
	(NAFC_lon,NAFC_lat),
	method='nearest'
	)

	#Check to see that index is correct
	for i,value in enumerate(NAFC_lon):
		if value < clim_lon[NAFC_ycoord[i],NAFC_xcoord[i]]:
			NAFC_xcoord[i] = NAFC_xcoord[i]-1
	for i,value in enumerate(NAFC_lat):
		if value < clim_lat[NAFC_ycoord[i],NAFC_xcoord[i]]:
			NAFC_ycoord[i] = NAFC_ycoord[i]-1

	#Determine which of the cast measurements are outliers
	outlier_mask = {}

	#Run through month by month
	for i in np.unique(ds['time.month'].values):

		#Import the climatology data
		path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/climatology/'
		#path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/climatology/'
		ds_clim = xr.open_dataset(path+'bootstrap_climatology_'+"%.2d" % (i)+'.nc')

		#Determine a monthly filter for the NAFC data
		month_filt = ds['time.month'].values == i
		outlier_mask[i] = {}

		#Cycle through temperature and salinity
		for ii in ['temperature','salinity']:

			#Create a variable to mask data
			outlier_mask[i][ii] = np.zeros((month_filt.sum(),ds_clim.level.size)).astype(bool)

			#Determine the climatology values at each cast location
			data_mean = ds_clim[ii+'_mean'].values[:,NAFC_ycoord[month_filt],NAFC_xcoord[month_filt]]
			data_ster = ds_clim[ii+'_ster'].values[:,NAFC_ycoord[month_filt],NAFC_xcoord[month_filt]]
			
			#Isolate the data
			data_NAFC = ds[ii][month_filt].values
			
			#Blanket 50 ster window for climatology
			place2 = \
			(data_NAFC > (data_mean.T + (data_ster.T*40)))+\
			(data_NAFC < (data_mean.T - (data_ster.T*40)))

			#Determine if there 25% of measurements are outliers
			nom_NAFC = data_NAFC.shape[1] - np.isnan(data_NAFC).sum(axis=1)
			nom_outliers = place2.sum(axis=1)
			np.seterr(divide='ignore', invalid='ignore')
			place2[nom_outliers/nom_NAFC > 0.25,:] = True
			outlier_mask[i][ii][:,:] = place2[:,:]

	#Create one universal filter for temperature and salinity
	outlier_combined = {}
	empties = {}
	max_outliers = {}

	#Cycle through temperature and salinity
	for ii in ['temperature','salinity']:

		#Create the empty mask
		outlier_combined[ii] = np.zeros((month_filt.size,ds_clim.level.size)).astype(bool)

		#Cycle through each month
		for i in np.unique(ds['time.month'].values):

			#Create the month filter, populate the outlier combined
			month_filt = ds['time.month'].values == i
			outlier_combined[ii][month_filt] = outlier_mask[i][ii]

		#Add back to the original combined dataset
		ds[ii][:,:] = ds[ii].values*(~outlier_combined[ii].astype(bool))

		#Determine where all the empties are
		empties[ii] = np.isnan(ds[ii]).sum(axis=1).values
		empties[ii] = empties[ii] == ds.level.size

		#Determine where the flagged cast for removal is 
		max_outliers[ii] = outlier_combined[ii].sum(axis=1) == ds.level.size


	#Determine where there are true empties (missing temperature and salinity)
	empties = empties['temperature']*empties['salinity']

	#Remove the cast with one or both variables exceeding the max number of outliers
	max_outliers = max_outliers['temperature']*max_outliers['salinity']

	#Combine the two
	cast_removal = empties+max_outliers

	#Change the temperature and salinity outliers to nans
	for ii in ['temperature','salinity']:
		place1 = ds[ii].values
		place1[(outlier_combined[ii])] = np.nan
		ds[ii][:,:] = place1

	#Remove the casts 
	ds = ds.sel(time = ~cast_removal)

	#Save the combined dataset
	#Flatten the arrays to ensure they save properly
	for i in ['instrument_ID','file_names','sounder_depth','instrument_type','station_ID','trip_ID','source']:
		ds[i] = ds[i].astype(str)
	#for i in ['trip_ID','comments','instrument_type','instrument_ID','file_names']:
	#	ds[i] = ds[i].astype(str)

	#Save the isolated netcdf files
	ds.to_netcdf('/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined_region_filtered/'+\
		str(year)+'.nc',
		encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	#ds.to_netcdf('/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/outliers/'+\
	#	str(year)+'.nc',
	#	encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	ds.close()
	print(year+' done.')


#Cycle through all the files and compress them
#Right now this needs to be done from nc_old after files have been moved there manually, not cool
for year in np.arange(1912,2022+1).astype(str)[:]:

	#Compress all files for memory purposes
	exp = 'nccopy -d 9 '+year+'.nc ../'+year+'.nc'
	os.system(exp)

	print(year+' done.')