import xarray as xr
import netCDF4 as nc
import numpy as np
import time as tt
import os
import glob
import pfile_tools as p
import pandas as pd
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt


#######################################################################################
'''

Remove duplicates from the .pfile sourced data.
Duplicates are defined as two or more casts which occur within the same hour and are
less than 0.01 degrees distance from one another.

The cast which contains more temperature recordings is kept.

'''
#######################################################################################


#Create a list of files made with the outliers removed
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/'
nc_list = [path+'outliers/'+i+'.nc' for i in np.arange(1912,2021+1).astype(str)]


#Cycle through each year
for i in np.sort(nc_list)[:]:

	#Load the .nc file, now with no empties
	pfile = xr.open_dataset(i)

	spatial_rounder = 2 #number of decimal points
	temporal_rounder = 'h' #hours

	#Determine the times where there are time duplicates
	time = pfile.time.values.astype('datetime64['+temporal_rounder+']')
	longitude = pfile.longitude.values.round(spatial_rounder)
	latitude = pfile.latitude.values.round(spatial_rounder)

	#Start to bring all the variables together
	spatiotemporal_comp = np.array([
		time.astype(str),
		latitude.astype(str),
		longitude.astype(str)
		])

	a = spatiotemporal_comp.T

	#Merge the time,lat,lon as one string, consider as a whole
	a = [row[0]+','+row[1]+','+row[2] for row in spatiotemporal_comp.T]
	unique_dates,unique_index,unique_count = np.unique(a,return_counts=True,return_index=True)

	#Mark the times as duplicates
	duplicates_flag = np.full(time.shape,False) 
	
	for ii in np.arange(unique_dates.size):

		#Determine if there is a count higher than 1 (duplicate)
		if unique_count[ii] > 1:

			#Isolate for the casts
			spots = np.where(np.array(a) == unique_dates[ii])[0]

			#Import the temperature
			temp = pfile.temperature[spots].values

			#Determine the number of recordings in each cast
			num_of_measurements = (~np.isnan(temp)).sum(axis=1)


			#Take the cast which has the higher number of measurements minus the number of outliers
			#In the case where the number of measurements is the same, take the first cast
			spots = spots[~(spots == spots[np.argmax(num_of_measurements)])]
			duplicates_flag[spots] = True


		#print(str(ii)+' out of '+str(duplicates_flag.size)+' done.')

	#Remove the duplicates from the pfile
	pfile = pfile.sel(time=~duplicates_flag)
	pfile.to_netcdf(path+'duplicates/'+i[-7:-3]+'.nc',mode='w')
	pfile.close()

	print(str(i[-7:-3])+' done.')


