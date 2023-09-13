import pfile_tools as p
import glob
import os
import warnings
import numpy as np
import xarray as xr
import netCDF4 as nc
import time as tt
import pandas as pd

'''

Port the NAFC .pfile files to an .nc file (1912-2022)
Base the structure of the .nc files from the yearly_netcdf_gen.py code

Charlie Bishop provides .pfiles (local ASCII file type) for each cast
New casts are provided each year from local monitoring programs (AZMP, AZOMP, etc.)
New casts are also pulled from the Marine Institute and other sources
Note. Starting ~2024, data will be provided in netcdf file types (structured according to CIOOS)
Script will need to be updated to manage the change in file type

Path, file names, and year variables are specific to source computer, change as necessary
Year variable can be changed to convert a specific year(s) pfiles

Note. pfile_tools package needed currently for this script
Be sure to add pfile_tools.py to your list of packages

'''


##########################################################################################################
##(1) - Create a netcdf file for each pfile

#Define the paths
path_pfiles = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/pfiles_1912-2022/'
path_lists = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/lists/'

#Define the years to cover
years = np.arange(1912,2022+1).astype(str)

#Run through and create lists
for year in years:

    #Use the find command for all years
    expr = 'find '+path_pfiles+year+'''/ -type f -iname '*.p'''+year+'''' > ./'''+year+'.list' 
    os.system(expr)
    expr = 'mv '+year+'.list '+path_lists
    os.system(expr)
    print(year)

#Create a path for where nc files will be saved
path_nc = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/pfiles_netcdf/'

#Cycle through each year and create individual nc files
for  year in years:

    #Import the list of files
    list_full = np.genfromtxt(path_lists+year+'.list', dtype=str)

    #Cycle through each of the files
    for file in list_full:

        #Determine file name and storage location
        outfile = path_nc+year+'/'+year+'_'+file[80:-6]+'.nc'
        warnings.simplefilter('ignore',category=RuntimeWarning)

        #Run Freds code
        try:
            p.pfiles_to_netcdf(file, outfile, zbin=1, zmax=5000)
        except:
            print('error occurred.')        





##########################################################################################################
##(2) - Merge individual netcdf files into yearly files

#Cycle through each year
for yearfile in np.arange(1912,2022+1).astype(str):

    #Merge the .nc files according to time and save into one .nc file
    #Create a str of the file names
    path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/pfiles_netcdf/'
    nc_list = np.array(glob.glob(str(path+yearfile)+'/*.nc'))

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
    density = np.zeros((nc_list.size,max_depth))

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
        density[i,:] = ds['sigma-t'].values[0,:]

        #Close the individual netcdf file
        ds.close()
        print(str(i)+' out of '+str(nc_list.size)+' done.')

    #Save the new variables as a netcdf array
    path_output = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/'
    nc_out = nc.Dataset(path_output+yearfile+'_test.nc','w')
 
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
    dens = nc_out.createVariable('density', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

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
    dens.units = 'Kg m-3'
    dens.long_name = 'Sigma-t'
    dens.standard_name = 'sigma_t'

    #Create a time filter
    filt = np.argsort(time_1D)

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
    dens[:,:] = density[filt,:]

    #Convert to time stamps
    time_stamps = [pd.Timestamp(i).to_pydatetime() for i in time_1D[filt]]
    times[:] = nc.date2num(time_stamps, units=times.units, calendar=times.calendar)
    levels[:] = np.arange(max_depth)

    #Save and close the .nc file
    nc_out.close()

    print(yearfile+' done.')



