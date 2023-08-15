'''

Port the BIO_Climate_Databases.RData file to an .nc file
Base the structure of the .nc files from the yearly_netcdf_gen.py code

This script converts both the Climate and BIO-OMO .RData files 
Note. Both source are discontinued, no updates are provided for each new year

Paths to .RData files are specific to sources computer, edit path variables as needed
For BIO-OMO, two files are provided (database_2008-2017.RData and database_2018.RData)
Change the file name variable in the second section to switch between 

'''

import xarray as xr
import netCDF4 as nc
import numpy as np
import time as tt
import os
import pfile_tools as p
import pandas as pd
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt
import pyreadr 

#Import one of the .RData files
file_name = 'BIO_Climate_Databases'
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CCAS20/'

#Load in the .RData file you want to convert
rfile = pyreadr.read_r(path+'R_Files/'+file_name+'.RData')



##########################################################################################################
##(1) - Pre-2010 Script, Used for BIO_Climate_Databases.RData

#Find out which years are covered by the file
years = np.arange(
	rfile['climate_stations'].CRUISE_DATE.values.min().astype('datetime64[Y]').astype(int)+1970,
	rfile['climate_stations'].CRUISE_DATE.values.max().astype('datetime64[Y]').astype(int)+1970+1)

#Cycle through each year
for year in years[:]:
	
	#Set up the .nc file
	nc_out = nc.Dataset(path+'NetCDF/'+file_name+'/'+str(year)+'_temporary.nc','w')
	
	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = '.RData test Conversion File' #Temporary title for the .nc file
	nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
	nc_out.source = 'https://github.com/jn533213/AZMP_Fred/blob/main/Trial2_RData.py'
	nc_out.references = 'Fill later'
	nc_out.description = 'A test by jonathan.coyne@dfo-mpo.gc.ca'
	nc_out.history = 'Created ' + tt.ctime(tt.time())

	#Create dimensions
	max_depth = 5000
	time = nc_out.createDimension('time', None) #use date2 for this
	level = nc_out.createDimension('level', max_depth) 

	#Create coordinate variables
	times = nc_out.createVariable('time', np.float64, ('time',))
	levels = nc_out.createVariable('level', np.int32, ('level',))

	#Create 1D variables
	cruise_id = nc_out.createVariable('Cruiseid', str, ('time'), zlib=True)
	latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
	longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
	cruisedate = nc_out.createVariable('cruisedate', str, ('time'), zlib=True)
	cruisetime = nc_out.createVariable('cruisetime', np.float32, ('time'), zlib=True)
	stnid = nc_out.createVariable('Stnid', np.float32, ('time'), zlib=True)
	datatype = nc_out.createVariable('datatype', str, ('time'), zlib=True)
	maximumdepth = nc_out.createVariable('maximumdepth', np.float32, ('time'), zlib=True)
	flag = nc_out.createVariable('flag', str, ('time'), zlib=True)
	samplecount = nc_out.createVariable('samplecount', np.float32, ('time'), zlib=True)
	update = nc_out.createVariable('update', str, ('time'), zlib=True)
	interpolated = nc_out.createVariable('interpolated', str, ('time'), zlib=True)

	#Create 2D variables
	#Each data variable is marked with a station ID
	#There are 13063 unique station IDs and 13089 station IDs in the rfile['stns']
	#Multiple data recordings exist for each unqiue station ID
	#Place each in the proper spot to create a 2D array and fill the empty values with nans
	temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	#pres = nc_out.createVariable('pressure', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	saln = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	sigm = nc_out.createVariable('sigmat', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	stni = nc_out.createVariable('stnid', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

	# Variable Attributes
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
	sigm.standard_name = "sigma_t"
	sigm.long_name = "Sigma-t"
	sigm.units = "Kg m-3"
	stni.standard_name = 'stnid' 
	stni.long_name = 'Station ID'

	#Create a filter for the specific year, will sort by time as well
	#Pre-2010
	filt = np.array([i.astype('datetime64[Y]').astype(int)+1970 for i in rfile['climate_stations'].CRUISE_DATE.values])
	filt = np.where(filt == year)[0]
	#Data is now sorted in terms of time 
	filt = filt[np.argsort(rfile['climate_stations'].CRUISE_DATE.values[filt])]

	#Fill in the 1D variables 
	cruise_id[:] = rfile['climate_stations'].CRUISE_ID.values[filt].astype(str)
	latitudes[:] = rfile['climate_stations'].LATITUDE.values[filt]
	longitudes[:] = rfile['climate_stations'].LONGITUDE.values[filt]
	cruisedate[:] = rfile['climate_stations'].CRUISE_DATE.values[filt].astype(str)
	cruisetime[:] = rfile['climate_stations'].CRUISE_TIME.values[filt]
	stnid[:] = rfile['climate_stations'].STN_ID.values[filt]
	datatype[:] = rfile['climate_stations'].DATATYPE.values[filt].astype(str)
	maximumdepth[:] = rfile['climate_stations'].MAXIMUM_DEPTH.values[filt]
	flag[:] = rfile['climate_stations'].FLAG.values[filt].astype(str)
	samplecount[:] = rfile['climate_stations'].SAMPLE_COUNT.values[filt]
	update[:] = rfile['climate_stations'].UP_DATE.values[filt].astype(str)
	interpolated[:] = rfile['climate_stations'].INTERPOLATED.values[filt].astype(str)


	#Re-format the 'data' variables into a 2D array
	temp2D = np.full((len(filt),max_depth),np.nan)
	saln2D = np.full((len(filt),max_depth),np.nan)
	sigm2D = np.full((len(filt),max_depth),np.nan)
	stni2D = np.full((len(filt),max_depth),np.nan)
	
	#Fill in the 2D arrays
	for i in rfile['climate_stations'].STN_ID[filt]:
		stn_spot = np.where(rfile['climate_stations'].STN_ID[filt] == i)[0][0] #Isolate the 'stns' spot
		dat_spot = np.where(rfile['climate_data'].STN_ID == i)[0] #Isolate all data instances of stnid

		dep_spot = rfile['climate_data'].PRESSURE.values[dat_spot].astype(int) #Isolate all relevant depths
		dat_spot = dat_spot[dep_spot < max_depth] #Isolate for depths shallower than max_depth 
		dep_spot = dep_spot[dep_spot < max_depth]    

		temp2D[stn_spot,dep_spot] = rfile['climate_data'].TEMPERATURE.values[dat_spot].astype(float)
		saln2D[stn_spot,dep_spot] = rfile['climate_data'].SALINITY.values[dat_spot].astype(float)
		sigm2D[stn_spot,dep_spot] = rfile['climate_data'].SIGMAT.values[dat_spot].astype(float)
		stni2D[stn_spot,dep_spot] = rfile['climate_data'].STN_ID.values[dat_spot].astype(float)

	#Filter 2D arrays for non-physical recordings and set to nan
	cutoff = {
	'temperature':	[-2,35],
	'salinity':		[0,45],
	}
	temp2D[temp2D <= cutoff['temperature'][0]] = np.nan
	temp2D[temp2D >= cutoff['temperature'][1]] = np.nan
	saln2D[saln2D <= cutoff['salinity'][0]] = np.nan
	saln2D[saln2D >= cutoff['salinity'][1]] = np.nan

	#Fill 2D structure
	temp[:,:] = temp2D
	saln[:,:] = saln2D
	sigm[:,:] = sigm2D
	stni[:,:] = stni2D

	#Convert to time stamps
	time_stamps = [pd.Timestamp(i).to_pydatetime() for i in rfile['climate_stations'].CRUISE_DATE.values[filt]]
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_depth)

	#Save and close the .nc file
	nc_out.close()
	print(str(year)+' done.')

#Cycle through each year and remove GTS-sources casts using the datatype variable
for year in years:

	#Open the dataset
	ds = xr.open_dataset(path+'NetCDF/'+file_name+'/'+str(year)+'_temporary.nc')

	#Add the file name variable
	ds['file_names'] = (('time'), np.tile(file_name,ds.time.size))
	#Remove the 'TE' and 'BA' datatype files (TESAC GTS and BATHY GTS)
	ds = ds.sel(time = ds.datatype != 'TE')
	ds = ds.sel(time = ds.datatype != 'BA')

	#Re-covert objects to strings
	for i in ['Cruiseid','cruisedate','datatype','flag','update','interpolated','file_names']:
		ds[i] = ds[i].astype(str)

	#Re-save the netcdf file
	if ds.time.size > 0:
		ds.to_netcdf(path+'NetCDF/'+file_name+'/'+str(year)+'.nc','w')
	ds.close()

	#Remove the temporary file
	expr = 'rm '+path+'NetCDF/'+file_name+'/'+str(year)+'_temporary.nc'
	os.system(expr)
	print(str(year)+' done.')





##########################################################################################################
##(2) - Post-2008 Script, Used for database_2008-2017.RData and database_2018.RData

#Import one of the .RData files
file_name = 'database_2018'
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/CCAS20/'

#Load in the .RData file you want to convert
rfile = pyreadr.read_r(path+'R_Files/'+file_name+'.RData')

#Find out which years are covered by the file
years = np.arange(
	np.datetime64(rfile['stns'].Date2.values[0]).astype('datetime64[Y]').astype(int)+1970,
	np.datetime64(rfile['stns'].Date2.values[-1]).astype('datetime64[Y]').astype(int)+1970+1)

#Cycle through each of the years
for year in years:
	
	#Create a .nc file
	nc_out = nc.Dataset(path+'NetCDF/'+file_name+'/'+str(year)+'_temporary.nc','w')

	#File information
	nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
	nc_out.title = '.RData test Conversion File' #Temporary title for the .nc file
	nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
	nc_out.source = 'https://github.com/jn533213/AZMP_Fred/blob/main/Trial2_RData.py'
	nc_out.references = 'Fill later'
	nc_out.description = 'A test by jonathan.coyne@dfo-mpo.gc.ca'
	nc_out.history = 'Created ' + tt.ctime(tt.time())

	#Create dimensions
	max_depth = 5000
	time = nc_out.createDimension('time', None) #use date2 for this
	level = nc_out.createDimension('level', max_depth) 

	#Create coordinate variables
	times = nc_out.createVariable('time', np.float64, ('time',))
	levels = nc_out.createVariable('level', np.int32, ('level',))

	#Create 1D variables
	cruise_id = nc_out.createVariable('Cruiseid', str, ('time'), zlib=True)
	latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
	longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
	cruisedate = nc_out.createVariable('cruisedate', str, ('time'), zlib=True)
	cruisetime = nc_out.createVariable('cruisetime', str, ('time'), zlib=True)
	stnid = nc_out.createVariable('Stnid', str, ('time'), zlib=True)
	datatype = nc_out.createVariable('datatype', str, ('time'), zlib=True)
	maximumdepth = nc_out.createVariable('maximumdepth', np.float32, ('time'), zlib=True)

	#Create 2D variables
	#Each data variable is marked with a station ID
	#There are 13063 unique station IDs and 13089 station IDs in the rfile['stns']
	#Multiple data recordings exist for each unqiue station ID
	#Place each in the proper spot to create a 2D array and fill the empty values with nans
	temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	#pres = nc_out.createVariable('pressure', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	saln = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	sigm = nc_out.createVariable('sigmat', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
	stni = nc_out.createVariable('stnid', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

	# Variable Attributes
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
	sigm.standard_name = "sigma_t"
	sigm.long_name = "Sigma-t"
	sigm.units = "Kg m-3"
	stni.standard_name = 'stnid' 
	stni.long_name = 'Station ID'

	#Create a filter for the specific year
	filt = np.array([np.datetime64(i).astype('datetime64[Y]').astype(int)+1970 for i in rfile['stns'].Date2.values])
	filt = np.where(filt == year)[0]
	#Data is now sorted in terms of time 
	filt_time = np.array([np.datetime64(i) for i in rfile['stns'].Date2.values])
	filt = filt[np.argsort(filt_time[filt])]

	#Post-2010
	cruise_id[:] = rfile['stns'].Cruiseid.values[filt]
	latitudes[:] = rfile['stns'].Latitude.values[filt]
	longitudes[:] = rfile['stns'].Longitude.values[filt]
	cruisedate[:] = rfile['stns'].Date2.values[filt]
	cruisetime[:] = rfile['stns'].Time.values[filt].astype(str)
	stnid[:] = rfile['stns'].Stnid.values[filt].astype(str)
	maximumdepth[:] = rfile['stns'].Pmax.values[filt]
	#Add datatype to it if 2018
	if file_name == 'database_2018':
		datatype[:] = rfile['stns']['Data Type'].values[filt].astype(str)

	#Re-format the 'data' variables into a 2D array
	temp2D = np.full((len(filt),max_depth),np.nan)
	saln2D = np.full((len(filt),max_depth),np.nan)
	sigm2D = np.full((len(filt),max_depth),np.nan)
	stni2D = np.full((len(filt),max_depth),np.nan)
	
	#Fill in the 2D arrays
	for i in rfile['stns'].Stnid[filt]:
		stn_spot = np.where(rfile['stns'].Stnid[filt] == i)[0][0] #Isolate the 'stns' spot
		dat_spot = np.where(rfile['data'].Stnid == i)[0] #Isolate all data instances of stnid

		dep_spot = rfile['data'].Pres.values[dat_spot].astype(int) #Isolate all relevant depths
		dat_spot = dat_spot[dep_spot < max_depth] #Isolate for depths shallower than max_depth 
		dep_spot = dep_spot[dep_spot < max_depth]    

		temp2D[stn_spot,dep_spot] = rfile['data'].Temp.values[dat_spot].astype(float)
		saln2D[stn_spot,dep_spot] = rfile['data'].Sal.values[dat_spot].astype(float)
		sigm2D[stn_spot,dep_spot] = rfile['data'].Sigm.values[dat_spot].astype(float)
		stni2D[stn_spot,dep_spot] = rfile['data'].Stnid.values[dat_spot].astype(float)

	#Filter 2D arrays for non-physical recordings and set to nan
	cutoff = {
	'temperature':	[-2,35],
	'salinity':		[0,45],
	}
	temp2D[temp2D <= cutoff['temperature'][0]] = np.nan
	temp2D[temp2D >= cutoff['temperature'][1]] = np.nan
	saln2D[saln2D <= cutoff['salinity'][0]] = np.nan
	saln2D[saln2D >= cutoff['salinity'][1]] = np.nan

	#Fill 2D structure
	temp[:,:] = temp2D
	saln[:,:] = saln2D
	sigm[:,:] = sigm2D
	stni[:,:] = stni2D

	#Convert to time stamps
	time_stamps = [pd.Timestamp(i).to_pydatetime() for i in filt_time[filt]]
	times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
	levels[:] = np.arange(max_depth)

	#Save and close the .nc file
	nc_out.close()
	print(str(year)+' done.')

#Cycle through each year and remove visually suspisious cruise id casts
for year in years:

	#Open the dataset
	ds = xr.open_dataset(path+'NetCDF/'+file_name+'/'+str(year)+'_temporary.nc')

	#Add the file name variable
	ds['file_names'] = (('time'), np.tile(file_name,ds.time.size))

	#Remove all casts with the following cruise_id starter
	#These casts appear to be buoy data
	cruise_id = ds.Cruiseid.values
	cruise_id_filter = np.full(cruise_id.size, False)
	for i,value in enumerate(cruise_id):

		#Check to see each of the starters
		if value.startswith('GoMoos'):
			cruise_id_filter[i] = True
		if value.startswith('OTN'):
			cruise_id_filter[i] = True
		if value.startswith('NERA'):
			cruise_id_filter[i] = True
		if value.startswith('79'):
			cruise_id_filter[i] = True
		if value.startswith('39'):
			cruise_id_filter[i] = True
		if value.startswith('12'):
			cruise_id_filter[i] = True

	#Filter out those casts
	ds = ds.sel(time = ~cruise_id_filter)
	ds = ds.sel(time = ds.datatype != 'TE') #For 2018
	ds = ds.sel(time = ds.datatype != 'BA') #For 2018

	#Remove all casts with 2 or less values
	ds = ds.sel(time = np.sum(~np.isnan(ds.temperature),axis=1) > 3)

	#Flatten the arrays to ensure they save properly
	for i in ['Cruiseid','cruisedate','cruisetime','Stnid','file_names','datatype']:
		ds[i] = ds[i].astype(str)

	#Re-save the netcdf file
	ds.to_netcdf(path+'NetCDF/'+file_name+'/'+str(year)+'.nc','w')
	ds.close()
	#Remove the temporary file
	expr = 'rm '+path+'NetCDF/'+file_name+'/'+str(year)+'_temporary.nc'
	os.system(expr)
	print(str(year)+' done.')






