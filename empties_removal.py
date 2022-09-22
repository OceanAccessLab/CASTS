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

Remove the empties from the .pfile sourced data
This will also include the pre-processing including removing un-physical measurements
and setting them as nans.

Casts which contain no measurements in any of the variables will be removed.

'''
#######################################################################################


#Create the .nc files which flag the empty casts
#Background info
source = 'netcdf_gen'
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/'

#Currently a problem with the naming of sigma-t, look into later

#First, create a list of .nc files in each path
#nc_list = glob.glob(path+'*.nc')
nc_list = [path+i+'.nc' for i in np.arange(1912,2021+1).astype(str)]

for i in np.sort(nc_list)[:]:

	##############################
	'''
	REMOVE EMPTIES 
	'''
	##############################

	#Set the nan values to missing
	year = i[-7:-3]

	##NEW VERSION USING XARRAY

	#Load in the dataset of interest
	ds = xr.open_dataset(i)

	#Run through each variable and determine the number of nans (missing) in each
	#Pre-define cut-offs for each variable that are non-realistic
	cutoff = {
	'temperature':	[-2,25],
	'salinity':		[0,45],
	'conductivity':	[-1e9,1e9],
	'oxygen':		[-1e9,1e9],
	'fluorescence':	[-1e9,1e9],
	'irradiance':	[-1e9,1e9],
	'ph':			[-1e9,1e9],
	}

	#Load in each of the variables one at a time
	missing = {}
	for ii in cutoff.keys():
		place1 = ds[ii].values

		#Perform the cutoffs listed above
		place1[place1 <= cutoff[ii][0]] = np.nan
		place1[place1 >= cutoff[ii][1]] = np.nan

		#Calculate the sum of missing values in each cast
		missing[ii] = np.isnan(place1).sum(axis=1)

	#Determine if a cast is completely empty
	missing = np.stack([missing[ii] for ii in cutoff.keys()])
	missing = missing.sum(axis=0) == missing.shape[0]*2000

	#Save the data now with the empty casts removed
	ds = ds.isel(time=~missing)
	ds.to_netcdf(path+'empties_duplicates/empties/'+year+'.nc')
	ds.close()

	print(year+', empty cast removal done.')

