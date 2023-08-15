import xarray as xr
import netCDF4 as nc
import numpy as np
import time as tt
import os
import csv
import subprocess
from datetime import datetime
import cftime
import pfile_tools as p
import pandas as pd
from scipy import stats,interpolate
import matplotlib
matplotlib.interactive(True)
import matplotlib.pyplot as plt
from urllib import request


'''

The purpose of this script is to convert the current structre of the Marine Institude
data in netcdf files to a consistent structure from other sources.

This data source is sparse and most likely incomplete
This source will not be updated in the future

Path inputs and outputs are specific to the source computer, change as needed

'''


##########################################################################################################
##(1) - Isolate and reformat casts according to year

#Define the path to the data
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Other/Freds_Files/'
path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Other/Freds_Files/data_processed'

#Define the maximum depth
max_depth = 5000

#Determine the file names
#folders = next(os.walk(path))[1]
folders = ['1995','1996','1998','1999','2000','2001','2002','2007','2008']

#Cycle through each folder
for folder in folders:

	#Determine the name of all files within each folder
	files = os.listdir(path+folder)

	#Cycle through each file
	for file in files:

		#Import the dataset
		ds = nc.Dataset(path+folder+'/'+file)

		#Determine if there are groups present
		if ds.groups == {}:

			#No groups present, skip
			continue

		#Get to the data, embedded in three folders...
		ds_sel = {}
		for i in ds.groups.keys():
			for ii in ds.groups[i].groups.keys():
				for iii in ds.groups[i].groups[ii].groups.keys():

					#Record the group of interest
					ds_sel[iii] = ds.groups[i].groups[ii].groups[iii]

		#Determine the number of casts
		noc = np.size(list(ds_sel.keys()))

		#Set up the 1D variables
		latitude_1D = np.zeros(noc)
		longitude_1D = np.zeros(noc)
		timestamp_1D = np.zeros(noc).astype('datetime64[m]')
		mission_id_1D = np.zeros(noc)

		#Set up the 2D variables
		dept_2D = np.full((noc,max_depth),np.nan)
		temp_2D = np.full((noc,max_depth),np.nan)
		cond_2D = np.full((noc,max_depth),np.nan)
		saln_2D = np.full((noc,max_depth),np.nan)
		sigm_2D = np.full((noc,max_depth),np.nan)

		#Cycle through each cast and record the data
		for i,value in enumerate(ds_sel):

			#Record the metadata
			latitude_1D[i] = ds_sel[value].latitude
			longitude_1D[i] = ds_sel[value].longitude
			timestamp_1D[i] = np.datetime64(cftime.num2date(ds_sel[value].timestamp, 'days since 1900-01-01'))
			mission_id_1D[i] = ds_sel[value].mission_id


			#Determine whether to use pressure or depth
			variables = list(ds_sel[value].variables)
			if np.isin(variables,['pres']).sum() == 0:
				call_level = 'depth'
			else:
				call_level = 'pres'

			#Record the 2D variables
			if np.isin(['temp','temp90'],variables).sum() > 0:
				call = np.array(['temp','temp90'])[np.isin(['temp','temp90'],variables)]
				temp_2D[i,:] = stats.binned_statistic(
					ds_sel[value].variables[call_level][:],
					ds_sel[value].variables[call[0]][:],
					'mean',
					bins=np.arange(max_depth+1)
					).statistic

			if np.isin(['saln','saln90'],variables).sum() > 0:
				call = np.array(['saln','saln90'])[np.isin(['saln','saln90'],variables)]
				saln_2D[i,:] = stats.binned_statistic(
					ds_sel[value].variables[call_level][:],
					ds_sel[value].variables[call[0]][:],
					'mean',
					bins=np.arange(max_depth+1)
					).statistic

			if np.isin(['cond','cond90'],variables).sum() > 0:
				call = np.array(['cond','cond90'])[np.isin(['cond','cond90'],variables)]
				cond_2D[i,:] = stats.binned_statistic(
					ds_sel[value].variables[call_level][:],
					ds_sel[value].variables[call[0]][:],
					'mean',
					bins=np.arange(max_depth+1)
					).statistic

			if np.isin(['depth'],variables).sum() > 0:
				call = np.array(['depth'])[np.isin(['depth'],variables)]
				dept_2D[i,:] = stats.binned_statistic(
					ds_sel[value].variables[call_level][:],
					ds_sel[value].variables[call[0]][:],
					'mean',
					bins=np.arange(max_depth+1)
					).statistic

			if np.isin(['sigmat'],variables).sum() > 0:
				call = np.array(['sigmat'])[np.isin(['sigmat'],variables)]
				sigm_2D[i,:] = stats.binned_statistic(
					ds_sel[value].variables[call_level][:],
					ds_sel[value].variables[call[0]][:],
					'mean',
					bins=np.arange(max_depth+1)
					).statistic

			#Implement cutoffs where applicable
			cutoff_temp = [-2,35]
			cutoff_saln = [0,45]

			temp_2D[temp_2D <= cutoff_temp[0]] = np.nan
			temp_2D[temp_2D > cutoff_temp[1]] = np.nan
			saln_2D[saln_2D <= cutoff_saln[0]] = np.nan
			saln_2D[saln_2D > cutoff_saln[1]] = np.nan

		#Save the file as a new formatted netcdf
		#Set up the .nc file
		nc_out = nc.Dataset(path_output+'/'+file,'w')
		
		#File information
		nc_out.Conventions = 'CF-1.6' #Ask Fred what this means
		nc_out.title = 'Marine Institute File-Specific NetCDF' #Temporary title for the .nc file
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
		mission_ids = nc_out.createVariable('mission_id', str, ('time'), zlib=True)
		file_names = nc_out.createVariable('file_names', str, ('time'), zlib=True)
		
		#Create 2D variables
		temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
		saln = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
		dens = nc_out.createVariable('density', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
		cond = nc_out.createVariable('conductivity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
		dept = nc_out.createVariable('depth', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

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
		dept.units = 'm'
		dept.long_name = 'Depth'
		dept.standard_name = 'depth'
		cond.units = 'NA'
		cond.long_name = 'Conductivity-90'
		cond.standard_name = 'conductivity_90'

		#Sort according to time
		filt = np.argsort(timestamp_1D)

		#Fill in the 1D variables 
		latitudes[:] = latitude_1D[filt]
		longitudes[:] = longitude_1D[filt]
		times[:] = timestamp_1D[filt]
		mission_ids[:] = mission_id_1D.astype(int).astype(str)[filt]
		file_names[:] = np.full(noc, file)

		#Fill 2D structure
		temp[:,:] = temp_2D[filt]
		saln[:,:] = saln_2D[filt]
		dens[:,:] = sigm_2D[filt]
		cond[:,:] = cond_2D[filt]
		dept[:,:] = dept_2D[filt]


		#Convert to time stamps
		time_stamps = [pd.Timestamp(i).to_pydatetime() for i in timestamp_1D[filt]]
		times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
		levels[:] = np.arange(max_depth)

		#Save and close the .nc file
		nc_out.close()
		print(file+' done.')


#Load in all of the files together
ds = xr.open_mfdataset(path_output+'/*.nc',combine='nested',concat_dim='time')

#Determine which years are available
years = ds['time.year'].values

#Cycle through each year
for year in np.unique(years):

	#Isolate the data for that year
	ds_isel = ds.sel(time = years == year)

	#Sort the data by time
	ds_isel = ds_isel.sortby('time')

	#Flatten the arrays to ensure they save properly
	for i in ['file_names','mission_id']:
		ds_isel[i] = ds_isel[i].astype(str)

	#Save and close the data
	ds_isel.to_netcdf('/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/Other/Freds_Files/netcdf_yearly/'+\
		str(year)+'.nc',
		encoding={'time':{'units': 'seconds since 1900-01-01 00:00:00', 'calendar': 'gregorian'}})
	ds_isel.close()
	print(str(year)+' done.')










