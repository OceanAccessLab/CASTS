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
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt


#######################################################################################
'''

Determine data records which are outliers (both spatially and temporally)
The spatial factor needs to consider the distance between data points

There will be some human-decided factors that need to be kept track of for this method

'''
#######################################################################################

#Background info
source = 'netcdf_gen'
variable_name = 'instrument_ID' #reference the header info

#Meshgrid lat/lon
#Define lat and lon on a constant grid across years
#According to the BIO_CLIMATE data paper latmin=35,latmax=80,lonmin=-100,lonmax=-42
dc = 1
x = np.arange(-100, -42, dc)
y = np.arange(35, 80, dc)  
lon, lat = np.meshgrid(x,y)

#Cycle through each year
years = np.arange(1912,2021+1)

for i in np.arange(years.shape[0])[1:]:

	#Load in the specific data
	path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/empties/'
	ds = xr.open_dataset(path+str(years[i])+'.nc')
	ds_time = np.array(ds['time.month'])
	ds_lons = np.array(ds.longitude)
	ds_lats = np.array(ds.latitude)
	ds_temp = np.array(ds.temperature)
	ds.close()

	#Create an empty array for all the binned measurements
	temp_3D = np.full((12,2000,y.shape[0]-1,x.shape[0]-1),np.nan)

	#Cycle through each of the months present
	for ii in np.unique(ds_time):

		#Cycle through each depth:
		for iii in np.arange(2000):

			#Determine if any measurements are present at this depth
			if np.isnan(ds_temp[:,iii][ds_time == ii]).all():
				None #No measurements present at this depth during this month

			else:
				#Bin the temperatures into the pre-defined latitudes and longitudes for specific months
				warnings.simplefilter("ignore", category=RuntimeWarning)
				place1 = stats.binned_statistic_2d(
					ds_lons[ds_time == ii],
					ds_lats[ds_time == ii],
					ds_temp[:,iii][ds_time == ii],
					np.nanmean,
					bins=[lon[0,:],lat[:,0]]
					).statistic

				#Record in the binned temperature
				temp_3D[ii-1,iii,:,:] = place1.T

	#Save the temperature for each month 
	#Cycle through each month to save
	for ii in np.arange(1,12+1):

		#Save the monthly climatology mean and stdv in a .nc format
		path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/climatology/'

		#Write a new .nc file
		nc_out = nc.Dataset(path+\
			'climatology_'+str(years[i])+'-'+"%.2d" % ii+'.nc','w')

		#File information
		nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
		nc_out.title = 'Gridded Mean temperature, NetCDF_Gen'  #Temporary title for the .nc file
		nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
		nc_out.source = 'https://github.com/jn533213/AZMP_Fred/blob/main/Trial2_RData.py'
		nc_out.references = 'Fill later'
		nc_out.description = 'A binned-temperature test by jonathan.coyne@dfo-mpo.gc.ca'
		nc_out.history = 'Created ' + tt.ctime(tt.time())

		#Create dimensions
		time = nc_out.createDimension('time', None) 
		level = nc_out.createDimension('level', 2000) 
		longitude = nc_out.createDimension('longitude', x.size-1) 
		latitude = nc_out.createDimension('latitude', y.size-1) 

		#Create coordinate variables
		times = nc_out.createVariable('time', np.float64, ('time',))
		levels = nc_out.createVariable('level', np.int32, ('level',))
		latitudes = nc_out.createVariable('latitude', np.float64, ('latitude',))
		longitudes = nc_out.createVariable('longitude', np.float64, ('longitude',))

		#Create the temperature variable
		temp = nc_out.createVariable(
			'temperature', 
			np.float64, 
			('time','level','latitude','longitude'), 
			zlib=True, fill_value=-9999)

		# Variable Attributes
		latitudes.units = 'degree_north'
		longitudes.units = 'degree_east'
		times.units = 'hours since 1900-01-01 00:00:00'
		times.calendar = 'gregorian'
		levels.units = 'dbar'
		levels.standard_name = "pressure"
		levels.valid_min = 0
		temp.units = 'Celsius'
		temp.long_name = "Water Temperature" # (may be use to label plots)
		temp.standard_name = "sea_water_temperature"

		#Save the temperature
		temp[0,:,:,:] = temp_3D[ii-1]

		#Convert to time stamps
		months = np.datetime64(str(years[i])+'-'+"%.2d" % ii)
		time_stamps = pd.Timestamp(months).to_pydatetime()

		times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
		levels[:] = np.arange(2000)
		latitudes[:] = y[:-1]
		longitudes[:] = x[:-1]

		nc_out.close()

	print(str(years[i])+' done.')


#Cycle through each month and determine the mean,stdv
#Best to reset the code here or you might run into memory issues
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/climatology/'
for i in np.arange(1,12+1):
	ds = xr.open_mfdataset(path+'climatology_*-'+"%.2d" % i+'.nc')

	#Determine the monthly average over the entire time period
	ds_mean = ds.mean('time')
	ds_stdv = ds.std('time')

	ds_mean.to_netcdf(path+'finished/climatology-mean-'+"%.2d" % i+'.nc')
	ds_mean.close()
	ds_stdv.to_netcdf(path+'finished/climatology-stdv-'+"%.2d" % i+'.nc')
	ds_stdv.close()
	ds.close()
	print(str(i)+' done.')


#######################################################################################
'''

This section of the script will flag outliers according to the previously calculated
climatology.

The criteria is based off of Gregory D.N., 2004 and UNESCO 1990. Note that the 
methods are not followed exactly as we do not go through the previous methods set 
by UNESCO 1990. 

Also, a sloppy landmask is made from HadISST in order to determine the distance of 
casts to land. Might be better to use a bathymetry data set but it's fine for now... 

How best to implement the outliers into the casts is still being discussed.
How many measurements in a cast should be flagged before the entire cast is removed?

'''
#######################################################################################


#Import the HadISST 
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/HadISST/'
ds = xr.open_dataset(path+'HadISST_sst.nc')

#Create the land mask - all nan values are land
land_mask = ds.sst[0].values
land_mask[~np.isnan(land_mask)] = 0.
land_mask[np.isnan(land_mask)] = 1.
land_mask[land_mask == 0.] = np.nan

#Determine all the land coordinates, Isolate for the North Atlantic
land_mask_lons,land_mask_lats = np.meshgrid(ds.longitude.values,ds.latitude.values) 
land_mask_lons = (land_mask_lons*land_mask)[0:70,70:160]
land_mask_lats = (land_mask_lats*land_mask)[0:70,70:160]

land_mask_coords = np.array([
land_mask_lons.flatten()[~np.isnan(land_mask_lons.flatten())],
land_mask_lats.flatten()[~np.isnan(land_mask_lats.flatten())]
])
ds.close()


#Create a mask for outliers for each year
#Background info
source = 'netcdf_gen'
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/'

#Declare path to climatology data
path_clim = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/climatology/'

#Cycle through each year
for year in np.arange(2002,2021+1):

	#Isolate the .nc data
	pfile = xr.open_dataset(path+'empties/'+str(year)+'.nc')
	time = np.array(pfile['time.month'])
	lons = np.array(pfile.longitude)
	lats = np.array(pfile.latitude)
	
	#Determine the distance to land for each cast
	distance_to_land = np.zeros(lons.size)
	R = 6373.0 #approx. radius of earth in km 

	for i in np.arange(lons.size):

		#The distance is calculated using the Haversine Formula, https://andrew.hedges.name/experiments/haversine/ 
		lat1 = math.radians(lats[i]) #Coordinates of casts
		lat2 = np.radians(land_mask_coords[1,:]) #Coordinates of land
		lon1 = math.radians(lons[i])
		lon2 = np.radians(land_mask_coords[0,:])

		dlon = lon2 - lon1
		dlat = lat2 - lat1

		a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
		c = 2*np.arctan2(np.sqrt(a), np.sqrt(1 - a))

		distance_to_land[i] = np.min(c*R)

	#Define a mask for the temperatures
	clim_mask = np.zeros((time.size,pfile.level.size),dtype=bool)

	#Cycle through each month
	for i in np.unique(time):

		month_filt = (pfile['time.month'] == i).values
		#Load in the temperature data for the month
		temp = pfile.temperature[month_filt].values

		#Import the corresponding climatology mean and stdv
		ds_mean = xr.open_dataset(path_clim+'finished/netcdf_gen-mean-'+"%.2d" % i+'_climatology.nc')
		ds_stdv = xr.open_dataset(path_clim+'finished/netcdf_gen-stdv-'+"%.2d" % i+'_climatology.nc')

		lon_mean,lat_mean = np.meshgrid(ds_mean.longitude.values, ds_mean.latitude.values)

		#Cycle through each level
		for ii in np.arange(pfile.level.size):
			#First, determine if the month is empty
			if temp.size == 0:
				None #Skip it

			else:
				#Second determine if the level is empty
				if np.isnan(temp[:,ii]).sum() == temp.shape[0]:
					clim_mask[month_filt][:,ii] = False

				else:
					#Interpolate the gridded structure to the lat and lon locations for the month
					depth_filt = np.isnan(temp[:,ii])
					mean = np.full(depth_filt.size,np.nan)
					stdv = np.full(depth_filt.size,np.nan)

					mean[~depth_filt] = interpolate.griddata(
						np.array((lon_mean.flatten(),lat_mean.flatten())).T,
						ds_mean.temperature[ii].values.flatten().T,
						(lons[month_filt][~depth_filt],lats[month_filt][~depth_filt])
						)
					stdv[~depth_filt] = interpolate.griddata(
						np.array((lon_mean.flatten(),lat_mean.flatten())).T,
						ds_stdv.temperature[ii].values.flatten().T,
						(lons[month_filt][~depth_filt],lats[month_filt][~depth_filt])
						)
					#Set 0 stdv to nan, issue with interpolate.griddata
					stdv[stdv <= 0] = np.nan

					#First isolate for casts within 1000km of land
					close_to_land = (distance_to_land < 1000)

					#Second, determine if the depth is less than 50m
					#If so, the climatology threshold is set at 5 times
					if ii <= 50:
						clim_mask[month_filt*(close_to_land),ii] = (np.abs(\
							temp[close_to_land[month_filt],ii]-\
							mean[close_to_land[month_filt]]) > (\
							stdv[close_to_land[month_filt]]*5))

						#For all others, set at 3 stdv
						#Determine if the absolute anomaly is greater than 3*stdv
						clim_mask[month_filt*(~close_to_land),ii] = (np.abs(\
							temp[~close_to_land[month_filt],ii]-\
							mean[~close_to_land[month_filt]]) > (\
							stdv[~close_to_land[month_filt]]*3))
					else:

						#If the depth is greater than 50m, always use 3*stdv
						clim_mask[month_filt,ii] = (np.abs(temp[:,ii] - mean ) > (stdv*3))

		ds_mean.close()
		ds_stdv.close()

	#Change the clim_mask into a nan (true), 1 (false) mask
	clim_mask = clim_mask.astype(float)
	clim_mask[clim_mask == 1.] = np.nan 
	clim_mask[clim_mask == 0.] = 1.

	#If there are over a certain number of outlier measurements in the cast, the whole cast is discarded
	pfile['temperature'] = pfile.temperature*clim_mask

	#If the number of outliers for a cast exceeds 10, remove the whole cast
	num_of_outliers = np.isnan(clim_mask).sum(axis=1)
	pfile = pfile.sel(time=num_of_outliers <= 10)

	#If the cast is now completely empty, remove
	empty_casts = np.isnan(pfile.temperature.values).sum(axis=1) == 2000
	pfile = pfile.sel(time=~empty_casts)

	#Compress the variables for file size management 
	comp = dict(zlib=True, complevel=5)
	encoding = {var: comp for var in pfile.data_vars}

	pfile.to_netcdf(path+'outliers/'+str(year)+'.nc',encoding=encoding)
	pfile.close()

	print(str(year)+' done.')






















#TESTING AREA
'''
plt.pcolor(
	ds_mean.longitude.values,
	ds_mean.latitude.values,
	ds_mean.temperature[0].values)
data = plt.scatter(
	lons[month_filt],
	lats[month_filt],
	c=stdv,s=30,edgecolors='Black',cmap='Reds')
plt.colorbar(data)


plt.figure()
for i in np.arange(1,13):
	ds_mean = xr.open_dataset(path+'finished/netcdf_gen-mean-'+"%.2d" % i+'_climatology.nc')
	ds_stdv = xr.open_dataset(path+'finished/netcdf_gen-stdv-'+"%.2d" % i+'_climatology.nc')
	plt.subplot(4,3,i)
	plt.pcolor(ds_mean.temperature[1],vmin=0,vmax=20)
	plt.title(str(i))
	plt.colorbar()

#Use this code, along with the Figures/temporary data to determine stdv of outliers per cast
outliers_stdv = []

for year in np.arange(pfiles['time.year'].min(),pfiles['time.year'].max()+1):
	num_of_outliers = np.load('num_of_outliers-'+str(year)+'.npy')

	plt.figure()
	plt.hist(num_of_outliers,bins=np.arange(1,30))
	plt.title(str(year)+' histogram')
	plt.savefig('histogram_'+str(year)+'.png')
	plt.close()

	num_of_outliers.astype(float)[num_of_outliers == 0] = np.nan
	outliers_stdv.append(np.nanstd(num_of_outliers))
	print(str(year)+' done.')
'''