import xarray as xr
import netCDF4 as nc
import numpy as np
from functools import reduce
import time as tt
import os
import csv
import pfile_tools as p
import pandas as pd
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt


'''

The purpose of this script is convert the CCGS Amundsen CTD data from the original 
.int file type to a uniform yearly formatted NetCDF file type.

'''


#Define the path of data
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Polar_Data_Catalogue/data_raw/'
path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Polar_Data_Catalogue/data_product/netcdf_group/'

#Set up conversion for the months (originally strings)
months = {
	'JAN':	'01','FEB':	'02',
	'MAR':	'03','APR':	'04',
	'MAY':	'05','JUN':	'06',
	'JUL':	'07','AUG':	'08',
	'SEP':	'09','OCT':	'10',
	'NOV':	'11','DEC':	'12',
}

#Determine names for all directories containing data
folders = next(os.walk(path))[1]

#Cycle through for each folder
for folder in np.sort(folders):


	#Determine the name of all files (including paths) within folder
	cast_files = []
	for paths, subdirs, files in os.walk(path+folder):
		for name in files:
			cast_files.append(os.path.join(paths, name))
	cast_files = np.array(cast_files)
	filt = np.zeros(cast_files.size).astype(bool)
	for i,value in enumerate(cast_files):
		if value[-3:] == 'int':
			filt[i] = True
		if value[-8:] == 'info.int':
			filt[i] = False
	cast_files = cast_files[filt]


	#Record the metadata for each cast
	cruise_number = np.zeros(cast_files.size).astype(str)
	cruise_name = np.zeros(cast_files.size).astype(str)
	file_name = np.zeros(cast_files.size).astype('<U64')
	time_1D = np.zeros(cast_files.size).astype('datetime64[s]')
	latitude = np.zeros(cast_files.size).astype(float)
	longitude = np.zeros(cast_files.size).astype(float)
	
	#Record the 2D variables of interest
	max_depth = 5000
	temperature = np.zeros((cast_files.size, max_depth)).astype(float)
	salinity = np.zeros((cast_files.size, max_depth)).astype(float)
	density = np.zeros((cast_files.size, max_depth)).astype(float)

	#Cycle through each cast file
	for i,cast_file in enumerate(cast_files):

		#Save the header
		with open(cast_file, encoding='ISO-8859-1') as line:
			reader = csv.reader(line, delimiter=':')
			cast_head = [row for row in reader][:]

		#Determine the header location
		#Determine the names of the pressure, temperature, salinity, density
		for ii,value in enumerate(cast_head):
			if value[0].split()[0] == 'Pres' or value[0].split()[0] == 'PRES' or value[0].split()[0] == 'Net':
				header_loc = ii

		#Pressure Titles	
		for ii,value in enumerate(cast_head):
			if value[-1] == ' Sea Pressure (sea surface - 0) [decibars]':
				pres_title = value[0][2:]
				break
			elif value[-1] == ' Net open pressure [db]':
				pres_title = value[0][2:]
				break
			pres_title = 'Pres'

		#Temperature Titles 
		for ii,value in enumerate(cast_head):
			if value[-1] == ' Temperature (1990 scale) [deg C]':
				temp_title = value[0][2:]
				break
			elif value[-1] ==  ' Mean Temperature (ITS-90) [degrees C]':
				temp_title = value[0][2:]
				break
			elif value[-1] ==  ' Temperature (ITS-90) [degrees C]':
				temp_title = value[0][2:]
				break				
			temp_title = 'Temp'

		#Salinity Titles
		for ii,value in enumerate(cast_head):
			if value[-1] == ' Practical Salinity [psu]':
				sal_title = value[0][2:]
				break
			elif value[-1] == ' Mean Salinity (PSS-78) [psu]':
				sal_title = value[0][2:]
				break
			elif value[-1] == ' Salinity (PSS-78) [psu]':
				sal_title = value[0][2:]
				break
			sal_title = 'Sal'

		#Density Titles 
		for ii,value in enumerate(cast_head):
			if value[-1] == ' Sigma-T (rho(s, t, 0)-1000) [kg/m^3]':
				den_title = value[0][2:]
				break
			den_title = 'Dens'

		#Import the cast data
		cast_data = pd.read_csv(
			cast_file,
			delim_whitespace=True, header=header_loc, encoding = 'ISO-8859-1')

		#Check to see if density is present
		if np.isin(den_title, cast_data.columns):
			density_present = True
		else:
			density_present = False

		#Record the 1D variables 
		for ii in cast_head:
			if ii[0] == '% Cruise_Number':
				cruise_number[i] = ii[1][1:]
			if ii[0] == '% Cruise_Name':
				cruise_name[i] = ii[1][1:]
			if ii[0] == '% Original_Filename':
				file_name[i] = folder+'__'+ii[1][1:]
			if ii[0] == '% Start_Date_Time [UTC]':
				try:
					#Sometimes months are string, sometimes number. Account for both
					time_1D[i] = np.datetime64(ii[1])
				except ValueError:
					second = ii [3]
					minute = ii[2]
					hour = ii[1][-2:]
					day = ii[1][1:3]
					month = months[ii[1][4:7].upper()]
					year =  ii[1][8:12]
					time_1D[i] = np.datetime64(year+'-'+month+'-'+day+'T'+hour+':'+minute+':'+second)
			if ii[0] == '% Initial_Latitude [deg]':
				latitude[i] = float(ii[1])
			if ii[0] == '% Initial_Longitude [deg]':
				longitude[i] = float(ii[1])

		#Record the 2D variables
		temperature[i,:] = stats.binned_statistic(
			cast_data[pres_title].values[1:].astype(float),
			cast_data[temp_title].values[1:].astype(float),
			'mean',
			bins = np.arange(max_depth+1)
			).statistic
		salinity[i,:] = stats.binned_statistic(
			cast_data[pres_title].values[1:].astype(float),
			cast_data[sal_title].values[1:].astype(float),
			'mean',
			bins = np.arange(max_depth+1)
			).statistic
		if density_present == True:
			density[i,:] = stats.binned_statistic(
				cast_data[pres_title].values[1:].astype(float),
				cast_data[den_title].values[1:].astype(float),
				'mean',
				bins = np.arange(max_depth+1)
				).statistic
		else:
			density[i,:] = np.nan

		#print(str(i)+' out of '+str(cast_files.size)+' done.')


	#Save in NetCDF format
	#Set up the .nc file
	nc_out = nc.Dataset(path_output+folder+'.nc','w')
	
	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'Polar Data Catalogue Yearly Netcdf File' #Temporary title for the .nc file
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
	cruise_numbers = nc_out.createVariable('cruise_number', str, ('time'), zlib=True)
	cruise_names = nc_out.createVariable('cruise_name', str, ('time'), zlib=True)
	file_names = nc_out.createVariable('file_names', str, ('time'), zlib=True)
	
	#Create 2D variables
	temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	saln = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	dens = nc_out.createVariable('density', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

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
	dens.units = 'Kg m-3'
	dens.long_name = 'Sigma-t'
	dens.standard_name = 'sigma_t'

	#Sort according to time
	filt = np.argsort(time_1D)

	#Fill in the 1D variables 
	latitudes[:] = latitude[filt]
	longitudes[:] = longitude[filt]
	times[:] = time_1D[filt]
	cruise_numbers[:] = cruise_number[filt]
	cruise_names[:] = cruise_name[filt]
	file_names[:] = file_name[filt]

	#Fill 2D structure
	temp[:,:] = temperature[filt]
	saln[:,:] = salinity[filt]
	dens[:,:] = density[filt]

	#Convert to time stamps
	time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_1D[filt]]
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_depth)

	#Save and close the .nc file
	nc_out.close()
	print(folder+' done.')


#Now cycle through all the NetCDF files and save by year
ds = {}
for folder in folders:
	ds[folder] = xr.open_dataset(path_output+folder+'.nc')

#Merge all of the sources together
ds_merged = xr.concat([ds[folder] for folder in folders],dim='time',combine_attrs='override')
ds_merged = ds_merged.sortby('time')

#Determine which years are covered
years = np.unique(ds_merged['time.year'])

#Cycle through each year
for year in years:

	#Slice the ds_merged by year
	ds_isolate = ds_merged.isel(time=ds_merged['time.year'] == year)
	ds_isolate.attrs = {
	'Conventions': 'CF-1.6',
	'title': 'Polar Data Catalogue Yearly NetCDF File',
	'institution': 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada',
	'source': 'https://github.com/OceanAccessLab/CASH',
	'references': 'Cyr et al., 2022',
	'description': 'Temperature and salinity cast records from the CASTS Dataset.',
	'history': 'Created ' + tt.ctime(tt.time())
	}

	#Save the merged dataset
	path_final = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Polar_Data_Catalogue/data_product/netcdf_yearly/'
	ds_isolate.to_netcdf(path_final+str(year)+'.nc',\
		encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	ds_isolate.close()
	print(str(year)+' done.')


