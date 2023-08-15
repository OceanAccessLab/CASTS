import xarray as xr
import netCDF4 as nc
import numpy as np
import time as tt
import warnings
import os
import math
import pfile_tools as p
import pandas as pd
from scipy import stats,interpolate
from scipy.signal import butter, lfilter, freqz
import scipy
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean


#######################################################################################
'''

Determine which outlier casts contain depths which are suspect according to a reference
dataset, in this case GEBCO.

This will be done using indexing rather than interpolation.
If the maximum recording of a cast exceeds 1.5 times the recorded bathymetry,
the cast is defined as suspect and removed.

'''
#######################################################################################

#Import the bathymetry data
ds_bath = xr.open_dataset('/gpfs/fs7/dfo/dpnm/joc000/Data/GEBCO_2022/gebco_2022_n80.0_s35.0_w-100.0_e-42.0.nc')

#Choose a resolution (by skipping) for the bathymetry
res = 5
bath_lon,bath_lat = np.meshgrid(ds_bath.lon[::res].values,ds_bath.lat[::res].values)
bath_elv = ds_bath.elevation[::res,::res].values.astype(float)
bath_elv[bath_elv >= 0] = np.nan

#Define the indices for the bathymetry
indx_lon,indx_lat = np.meshgrid(np.arange(bath_lon.shape[1]),np.arange(bath_lat.shape[0]))


#Define the path and years
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined_region/'
path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Products/combined_depth_filtered/'
years = np.arange(1912,2022+1).astype(str)

#Cycle through each year
for year in years:

	#Open up the relevant data set
	ds = xr.open_dataset(path+year+'.nc')

	#Record all the latitudes and longitude
	ds_lon = ds.longitude.values
	ds_lat = ds.latitude.values

	#Determine the maximum depth of the casts
	#Load in the temperature and salinity
	temp = ds.temperature.values
	saln = ds.salinity.values

	#Determine the depth of each cast (deepest measurement)
	depth_max_temp = np.zeros(ds.time.size)
	depth_max_saln = np.zeros(ds.time.size)
	level = ds.level.values

	#Cycle through each cast
	for i in np.arange(ds.time.size):

		#Determine where the last measurement is
		if np.isnan(temp[i]).all() == False:
			depth_max_temp[i] = level[np.where(np.isnan(temp[i,:]) == False)[0][-1]]
		if np.isnan(saln[i]).all() == False:
			depth_max_saln[i] = level[np.where(np.isnan(saln[i,:]) == False)[0][-1]]

	#Determine where the temperature or salinity is the maximum
	depth_max = np.zeros(ds.time.size)
	place1 = depth_max_temp > depth_max_saln
	depth_max[place1] = depth_max_temp[place1]
	depth_max[~place1] = depth_max_saln[~place1]

	#Determine the index in bathymetry for each of the casts
	ds_xcoord = interpolate.griddata(
	np.array((bath_lon.flatten(),bath_lat.flatten())).T,
	indx_lon.flatten().T,
	(ds_lon,ds_lat),
	method='nearest'
	)
	ds_ycoord = interpolate.griddata(
	np.array((bath_lon.flatten(),bath_lat.flatten())).T,
	indx_lat.flatten().T,
	(ds_lon,ds_lat),
	method='nearest'
	)

	#Determine the GEBCO bathymetry at each cast point
	depth_GEBCO = np.abs(bath_elv[ds_ycoord,ds_xcoord])

	#Determine whether casts are exceeding the depth threshold
	#If the cast recorded depth is greater than 200m
	#If the cast depth is greater than 1.5 times the GEBCO depth
	#Then remove cast
	depth_filter = np.zeros(depth_max.size).astype(bool)
	for i in np.arange(depth_filter.size):
		if depth_max[i] > 200:
			if (depth_max[i] - depth_GEBCO[i]) > 300:
				depth_filter[i] = True


	#Save the casts with the depth-filtered casts removed
	ds = ds.sel(time = ~depth_filter)

	#Flatten the arrays to ensure they save properly
	for i in ['trip_ID','source','instrument_ID','instrument_type','file_names','station_ID','sounder_depth']:
		ds[i] = ds[i].astype(str)
	ds.to_netcdf(path_output+year+'.nc')
	ds.close()

	print(year+' done.')




