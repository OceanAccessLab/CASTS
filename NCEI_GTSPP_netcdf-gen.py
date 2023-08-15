import xarray as xr
import netCDF4 as nc
import numpy as np
import time as tt
import os
import csv
import subprocess
import pfile_tools as p
import pandas as pd
from scipy import stats,interpolate
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt
from urllib import request


'''

The purpose of this script is to convert the current structre of the NCEI_GTSPP files 
(.meds) into yearly netcdf output for combination with the other data sets late.

This script includes a programmed pull of data from the source website, this can take time
The path inputs and outputs are specific to the source computer, change as needed

This script also requires a GEBCO bathymetry netcdf file in order to perform step 3
Please download the appropriate file from the link below
https://download.gebco.net/
The following boundaries were used: 35-80 deg N, 42-100 deg W

'''


##########################################################################################################
##(1) - Pull the relevant casts from NOAA


#Pull the relevant files from HTTP addresses 
#This is meant to be run locally on the gpsc

#Define the result path, years, months
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NCEI_GTSPP/data_raw/'
years = np.arange(1985,2022+1).astype(str)
months = np.array(['01','02','03','04','05','06','07','08','09','10','11','12'])

#Define the web address
address = 'https://www.ncei.noaa.gov/data/oceans/gtspp/bestcopy/netcdf/'

#Cycle through each year
for year in years:

	#Create a yearly folder
	exp = 'mkdir '+path+year
	os.system(exp)

	#Cycle through each month
	for month in months:

		#Create a monthly folder
		exp = 'mkdir '+path+year+'/'+month
		os.system(exp)

		#Import the file, not always present
		try:
			request.urlretrieve(address+'gtspp4_at'+year+month+'.tgz',
				path+year+'/'+'gtspp4_at'+year+month+'.tgz')		
		except:
			continue

		#Open the file
		exp = 'tar -xvzf '+path+year+'/gtspp4_at'+year+month+'.tgz -C '+path+year+'/'
		os.system(exp)

		#Move the relevant files to the monthly folder
		#data_types = ['bo','cd','ct','cu','bf','to','xc','mb','dt','xb']
		#for i in data_types:
		exp = 'echo '+path+year+'/atlantic/'+year+'/'+month+'/*.nc | xargs mv -t '+path+year+'/'+month+'/ --'
		os.system(exp)

		#Delete the atlantic folder, move tgz file
		exp = 'rm -r '+path+year+'/atlantic'
		os.system(exp)
		exp = 'mv '+path+year+'/gtspp4_at'+year+month+'.tgz '+path+'tgz_files/'
		os.system(exp)

		print(year+'-'+month+' done.')


##########################################################################################################
##(2) - Compile like months together into netcdf files

#Cycle through each month and compile the months together into on .nc file
#Define spatial domain limits
lon_min = -100
lon_max = -42
lat_min = 35
lat_max = 80

#Cycle through each year
for year in years[27:]:

	#Cycle through each month
	for month in months:

		#Create a list of the netcdf files
		nc_list = np.array(os.listdir(path+year+'/'+month))
		max_depth = 5000

		#Determine which files are GTS (TE and BA)
		GTS_filt = np.zeros(nc_list.size).astype(bool)
		if nc_list.size > 0:
			for i,name in enumerate(nc_list):
				if name[-9:-7] == 'te' or name[-9:-7] == 'ba':
					GTS_filt[i] = True
			nc_list = nc_list[~GTS_filt]

		#Determine if any casts are present for the month
		if nc_list.size > 0:

			#Create empty 1D variables
			time_1D = np.zeros(nc_list.size).astype('datetime64[s]')
			latitude_1D = np.zeros(nc_list.size)
			longitude_1D = np.zeros(nc_list.size)
			data_type_1D = np.zeros(nc_list.size).astype(str)
			file_name_1D = np.zeros(nc_list.size).astype(str)
			cruise_id_1D = np.zeros(nc_list.size).astype(str)
			source_id_1D = np.zeros(nc_list.size).astype(str)

			#Create empty 2D variables
			temp_2D = np.zeros((nc_list.size, max_depth))
			saln_2D = np.zeros((nc_list.size, max_depth))

			#Cycle through each cast
			for i,cast in enumerate(nc_list):

				#Open the cast dataset
				ds = xr.open_dataset(path+year+'/'+month+'/'+cast)

				#Populate the 1D variables
				time_1D[i] = ds.time.values[0]
				latitude_1D[i] = ds.latitude.values[0]
				longitude_1D[i] = ds.longitude.values[0]
				data_type_1D[i] = ds.data_type.values.astype(str)
				cruise_id_1D[i] = ds.cruise_id.values.astype(str)
				source_id_1D[i] = ds.source_id.values.astype(str)
				file_name_1D[i] = cast

				#Populate the 2D variables
				#Check to see if salinity is available
				if np.isin('temperature',list(ds.variables)) == True:
					if np.isin(ds.temperature_quality_flag.values, [1]).sum() > 0:
						temp_2D[i,:] = stats.binned_statistic(
							ds.z.values[np.isin(ds.temperature_quality_flag.values, [1])],
							ds.temperature.values.flatten()[np.isin(ds.temperature_quality_flag.values, [1])],
							'mean',
							bins = np.arange(max_depth+1)
							).statistic
					else:
						temp_2D[i,:] = np.nan
				else:
					temp_2D[i,:] = np.nan

				#Check to see if salinity is available
				if np.isin('salinity',list(ds.variables)) == True:
					if np.isin(ds.salinity_quality_flag.values, [1]).sum() > 0:
						saln_2D[i,:] = stats.binned_statistic(
							ds.z.values[np.isin(ds.salinity_quality_flag.values, [1])],
							ds.salinity.values.flatten()[np.isin(ds.salinity_quality_flag.values, [1])],
							'mean',
							bins = np.arange(max_depth+1)
							).statistic
					else:
						saln_2D[i,:] = np.nan
				else:
					saln_2D[i,:] = np.nan

				#Close the dataset
				ds.close()
				print(str(i)+' out of '+str(nc_list.size)+' done. '+year+'-'+month+'.')


			#Create a nc file for all the casts
			#Set up the individual .nc file
			nc_out = nc.Dataset(path+year+'/'+year+'_'+month+'_temporary.nc','w')

			#File information
			nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
			nc_out.title = 'NCEI_GTSPP, NetCDF File' #Temporary title for the .nc file
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
			data_types = nc_out.createVariable('data_type', str, ('time'), zlib=True)
			cruise_IDs = nc_out.createVariable('cruise_ID', str, ('time'), zlib=True)
			file_names = nc_out.createVariable('file_name', str, ('time'), zlib=True)
			source_IDs = nc_out.createVariable('source_ID', str, ('time'), zlib=True)


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

			#Sort the variables according to time
			filt = np.argsort(time_1D)

			#Fill in the 1D variables
			latitudes[:] = latitude_1D[filt]
			longitudes[:] = longitude_1D[filt]
			times[:] = time_1D[filt]
			cruise_IDs[:]  = cruise_id_1D[filt]
			data_types[:] = data_type_1D[filt]
			file_names[:] = file_name_1D[filt]
			source_IDs[:] = source_id_1D[filt]

			#Fill 2D structure
			temp[:,:] = temp_2D[filt,:]
			saln[:,:] = saln_2D[filt,:]

			#Convert to time stamps
			time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_1D[filt]]
			times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
			levels[:] = np.arange(max_depth)

			#Save and close the .nc file
			nc_out.close()

			#Open the file using xarray to isolate for latitude and longitude range
			ds = xr.open_dataset(path+year+'/'+year+'_'+month+'_temporary.nc')

			#Isolate for the lats and lons
			ds = ds.sel(time=(ds.longitude > lon_min)*(ds.longitude < lon_max))
			ds = ds.sel(time=(ds.latitude > lat_min)*(ds.latitude < lat_max))

			#Re-covert objects to strings
			for i in ['data_type','file_name','cruise_ID','source_ID']:
				ds[i] = ds[i].astype(str)

			#Save the new netcdf
			if ds.time.size > 0:
				ds.to_netcdf(path+year+'/'+year+'_'+month+'_combined.nc','w')
				ds.close()

			#Remove the temporary file
			exp = 'rm '+path+year+'/'+year+'_'+month+'_temporary.nc'
			os.system(exp)

	#Now that the monthly files are made, convert into one yearly file
	path_output = path+'data_processed/'

	#Open the files together
	ds = xr.open_mfdataset(path+year+'/'+year+'_*_combined.nc')
	
	#Remove all GTS-sourced casts
	ds = ds.sel(time = ds.data_type != 'TE')
	ds = ds.sel(time = ds.data_type != 'BA')
	
	#Re-covert objects to strings
	for i in ['data_type','file_name','cruise_ID','source_ID']:
		ds[i] = ds[i].astype(str)

	#Save the new netcdf
	ds.to_netcdf(path_output+year+'.nc','w')
	ds.close()


##########################################################################################################
##(3) - Depth check to ensure that casts measurements are reasonable


#There has been issues with casts from this source going past reasonable depths
#Do a double check here, determine the bathymetry and check with the lowest measurement of each cast
#This takes a bit of time since each point needs to be interpolated

#Load in the bathymetry data
ds_bath = xr.open_dataset('/gpfs/fs7/dfo/dpnm/joc000/Data/GEBCO_2022/gebco_2022_n80.0_s35.0_w-100.0_e-42.0.nc')

#Interpolate the bathymetry to each of the cast locations
#Choose a resolution (by skipping) for the bathymetry
res = 10
lons,lats = np.meshgrid(ds_bath.lon.values[::res],ds_bath.lat.values[::res])
points = np.array([lons.flatten(),lats.flatten()]).T
values = ds_bath.elevation[::res,::res].values


#Load in the data, one year at a time
for year in years[5:]:
	ds = xr.open_dataset(path+'data_processed/'+year+'.nc')

	#Load in the latitude and longitude
	ds_lat = ds.latitude.values
	ds_lon = ds.longitude.values

	#Load in the temperature and salinity
	temp = ds.temperature.values
	saln = ds.salinity.values

	#Determine the depth of each cast (deepest measurement)
	depth_max_temp = np.zeros(ds.time.size)
	depth_max_saln = np.zeros(ds.time.size)

	#Cycle through each cast
	for i in np.arange(ds.time.size):

		#Determine where the last measurement is
		if np.isnan(ds.temperature[i,:].values).all() == False:
			depth_max_temp[i] = np.where(np.isnan(temp[i,:]) == False)[0][-1]
		if np.isnan(ds.salinity[i,:].values).all() == False:
			depth_max_saln[i] = np.where(np.isnan(saln[i,:]) == False)[0][-1]

	#Determine where the temperature or salinity is the maximum
	depth_max = np.zeros(ds.time.size)
	place1 = depth_max_temp > depth_max_saln
	depth_max[place1] = depth_max_temp[place1]
	depth_max[~place1] = depth_max_saln[~place1]

	#Determine the depth at each point
	bath_depth = interpolate.griddata(points, values.flatten(), (ds_lon,ds_lat),method='linear')

	#Flag the casts that are too deep and remove
	bottom_flag = (depth_max*-1) < bath_depth
	ds = ds.sel(time = ~bottom_flag)

	#Re-covert objects to strings
	for i in ['data_type','file_name','cruise_ID','source_ID']:
		ds[i] = ds[i].astype(str)

	#Save the new file
	ds.to_netcdf(path+'data_processed/bottom_flagged/'+year+'.nc')
	ds.close()
	print(year+' done.')






































