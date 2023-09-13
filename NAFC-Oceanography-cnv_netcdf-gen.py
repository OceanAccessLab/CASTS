import glob
import os
import sys
import warnings
import numpy as np
from scipy import stats,interpolate
import xarray as xr
import netCDF4 as nc
import time as tt
import pandas as pd
sys.path.append('/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Other/Charlies_Scripts')
import cnv_tk
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt

'''

Port the NAFC .pcnv files to an .nc file (2023-)
Base the structure of the .nc files from the yearly_netcdf_gen.py code

Charlie Bishop provides .pcnv (temporary replacement for .pfile) for each cast
New casts are provided each year from local monitoring programs (AZMP, AZOMP, etc.)
New casts are also pulled from the Marine Institute and other sources
Note. Starting ~2024, data will be provided in netcdf file types (structured according to CIOOS)
Script will need to be updated to manage the change in file type

Path, file names, and year variables are specific to source computer, change as necessary
Year variable can be changed to convert a specific year(s) pfiles

'''


##########################################################################################################
##(1) - Create a netcdf file for each pfile

#Define the paths
path_pcnv = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/pfiles_1912-2022/'
path_nc = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/pfiles_netcdf/'

#Define the years to cover
years = np.arange(2023,2024).astype(str)

#Cycle through each year and create individual nc files
for  year in years:

	#Determine all the files present
	file_names = os.listdir(path_pcnv+year)

	#Cycle through each of the files, save as netcdf
	for file_name in file_names:

		#Create the cast object
		cast = cnv_tk.Cast(path_pcnv+year+'/'+file_name)
		cnv_tk.cnv_meta(cast, path_pcnv+year+'/'+file_name)
		#Create a dataframe
		df = cnv_tk.cnv_to_dataframe(cast)
		df = cnv_tk.df_press_depth(cast)
		df = cnv_tk.StandardizedDF(cast, df)

		#Bin average the variables in 1dbar bins
		max_depth = 5000
		variables = ['Temperature','Salinity']
		variables_binned = {}
		for variable in variables:
			if variable in df.columns:
				variables_binned[variable] = stats.binned_statistic(
					df['Pressure'].values.astype(float),
					df[variable].values.astype(float),
					'mean',
					bins=np.arange(max_depth+1)
					).statistic
			else:
				variables_binned[variable] = np.full(max_depth,np.nan)



		#Save as a netcdf file
		nc_out = nc.Dataset(path_nc+year+'/'+file_name.split('.')[0]+'.nc','w')

		#File information
		nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
		nc_out.title = 'NAFC-Oceanography, NetCDF File' #Temporary title for the .nc file
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
		trip_IDs = nc_out.createVariable('trip_ID', str, ('time'), zlib=True)
		comments = nc_out.createVariable('comments', str, ('time'), zlib=True)
		instrument_types = nc_out.createVariable('instrument_type', 'U13', ('time'), zlib=True)
		instrument_IDs = nc_out.createVariable('instrument_ID', str, ('time'), zlib=True)
		sounder_depths = nc_out.createVariable('sounder_depth', str, ('time'), zlib=True)

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

		#Fill in the 1D variables
		latitudes[:] = cast.Latitude
		longitudes[:] = cast.Longitude
		times[:] = pd.Timestamp(cast.CastDatetime).to_datetime64()
		instrument_types[:] = np.array([cast.castType])
		instrument_IDs[:] = np.array([cast.InstrumentName])
		sounder_depths[:] = np.array([cast.SounderDepth])
		trip_IDs[:] = np.array([cast.triptag+'_'+cast.trip])
		comments[:] = np.array([cast.comment])

		#Fill in the 2D variables
		temp[:,:] = variables_binned['Temperature'][None,:]
		saln[:,:] = variables_binned['Salinity'][None,:]

		#Convert to time stamps
		time_stamps = pd.Timestamp(cast.CastDatetime).to_pydatetime()
		times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
		levels[:] = np.arange(max_depth)

		#Save and close the .nc file
		nc_out.close()
		print(file_name+' done.')



##########################################################################################################
##(2) - Merge individual netcdf files into yearly files

#Cycle through each year
for year in np.arange(2023,2024).astype(str):

	#Merge the .nc files according to time and save into one .nc file
	#Create a str of the file names
	path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/pfiles_netcdf/'
	nc_list = np.array(glob.glob(str(path+year)+'/*.nc'))

	#Open the files and save the necessary data
	#Set up the 1D arrays
	latitude = np.zeros(nc_list.size)
	longitude = np.zeros(nc_list.size)
	time_1D = np.zeros(nc_list.size).astype('datetime64[s]')
	trip_id = np.zeros(nc_list.size).astype(str)
	comment = np.zeros(nc_list.size).astype(str)
	instrument_type = np.zeros(nc_list.size).astype(str)
	instrument_id = np.zeros(nc_list.size).astype(str)
	file_name = np.zeros(nc_list.size).astype(str)

	#Set up the 2D arrays
	max_depth = 5000
	temperature = np.zeros((nc_list.size,max_depth))
	salinity = np.zeros((nc_list.size,max_depth))

	#Cycle through each nc file
	for i,file in enumerate(nc_list):

		#Open the file
		ds = xr.open_dataset(file)

		#Populate the 1D arrays
		latitude[i] = ds.latitude.values[0]
		longitude[i] = ds.longitude.values[0]
		time_1D[i] = ds.time.values[0]
		trip_id[i] = ds.trip_ID.values[0]
		comment[i] = ds.comments.values[0]
		instrument_type[i] = ds.instrument_type.values[0]
		instrument_id[i] = ds.instrument_ID.values[0]
		file_name[i] = file.split('/')[-1]

		#Populate the 2D arrays
		temperature[i,:] = ds.temperature.values[0,:]
		salinity[i,:] = ds.salinity.values[0,:]

		#Close the individual netcdf file
		ds.close()
		print(str(i)+' out of '+str(nc_list.size)+' done.')

	#Save the new variables as a netcdf array
	path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/'
	nc_out = nc.Dataset(path_output+year+'_test.nc','w')

	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = 'NAFC Yearly Netcdf File' #Temporary title for the .nc file
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
	trip_ids = nc_out.createVariable('trip_ID',str,('time'),zlib=True)
	comments = nc_out.createVariable('comments',str,('time'),zlib=True)
	instrument_types = nc_out.createVariable('instrument_type',str,('time'),zlib=True)
	instrument_ids = nc_out.createVariable('instrument_ID',str,('time'),zlib=True)
	file_names = nc_out.createVariable('file_names',str,('time'),zlib=True)    

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
	#levels.valid_min = 0
	temp.units = 'Celsius'
	temp.long_name = "Water Temperature" # (may be use to label plots)
	temp.standard_name = "sea_water_temperature"
	saln.long_name = "Practical Salinity"
	saln.standard_name = "sea_water_salinity"
	saln.units = "1"
	#saln.valid_min = 0

	#Create a time filter
	filt = np.argsort(time_1D)
	filt = filt[time_1D >= np.datetime64(year+'-01-01')]

	#Fill in the 1D variables 
	latitudes[:] = latitude[filt]
	longitudes[:] = longitude[filt]
	times[:] = time_1D[filt]
	trip_ids[:] = trip_id[filt]
	comments[:] = comment[filt]
	instrument_types[:] = instrument_type[filt]
	instrument_ids[:] = instrument_id[filt]
	file_names[:] = file_name[filt]

	#Fill 2D structure
	temp[:,:] = temperature[filt,:]
	saln[:,:] = salinity[filt,:]

	#Convert to time stamps
	time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_1D[filt]]
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_depth)

	#Save and close the .nc file
	nc_out.close()

	print(year+' done.')
