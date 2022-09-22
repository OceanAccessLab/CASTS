"""A script to loop on pfile lists and generate multiple (yearly) netCDF files


I need to make it a function!


** Note; works also for Viking buoy Dfiles.
Just generate the list with:
ls viking_pfiles/*.D20 > viking2020.list

----------
Frederic.Cyr@dfo-mpo.gc.ca, October 2017

"""

# NOTE FOR J COYNE:
# 1. In a terminal in HP500/netcdf_gen' do: for i in `seq 2011 2021`; do ls pfiles/*.p$i > ./$i.list; done
# (will take time! and may cause memory errors... see alternate method commented below)
# 2. In the same folder, open python and do: %run yearly_netcdf_gen.py (this file!)
# ** you need to make sure to retrive "ptools" here: github.com/AZMP-NL/python-toolbox/tree/master/ptools

import pfile_tools as p
import glob
import os
import warnings
import numpy as np
import xarray as xr
import netCDF4 as nc

#I expect to see warnings from the following loop, Jon C.
#RuntimeWarning, mean of empty slice is usually repeating throughout

#Deal with filenames containing spaces, some from Bishop do
#Replace ' ' with '_', cd into the pfiles folder!
#expr = 'for f in *\ *; do mv "$f" "${f// /_}"; done'
#os.system(expr)

#Create lists for each of the years of each of the pfiles
expr = 'for i in `seq 1912 2010`; do ls pfiles/*.p$i > ./$i.list; done'
os.system(expr)

lists = glob.glob('*.list')
lists = np.sort(lists)

os.system('echo " --------- New run ------ " >> .netcdfgen_log.txt')


#Use the cdo module to combine the pfile .nc files, testing phase right now

#Cycle through the lists
for yearfile in lists:

    #Create a .nc file for each component of the list
    list_full = np.genfromtxt(yearfile,dtype=str,replace_space='/ ')

    #Create a directory for all the individual .pfiles
    expr = 'mkdir '+str(yearfile[:4])
    os.system(expr)

    #Create the individual .nc files for each .pfile
    for i in list_full:

        outfile = str(yearfile[:4])+'/'+str(yearfile[:4])+'_'+i[7:-6]+'.nc'
        warnings.simplefilter('ignore',category=RuntimeWarning)

        try:
            p.pfiles_to_netcdf(i, outfile, zbin=1, zmax=2000)
            '''
            #Delete the file if the data is blank
            if xr.open_dataset(outfile).time.values == :
                expr = 'rm '+outfile
                os.system(expr)
            '''

        except:
            print('error occurred.')
            '''
            #Create a text document to catologue the errors
            expr = 'echo error occurred, '+i+' >> .netcdfgen_log.txt'
            os.system(expr)
            '''


    #Merge the .nc according to time and save under one .nc
    #Account for years with larger list sizes
    if list_full.size > 1000:

        #Create a str of the file names
        nc_list = glob.glob(str(yearfile[:4])+'/*.nc')

        for i in np.arange(0,list_full.size,1000):

            expr = 'cdo mergetime '+\
            ' '.join(nc_list[i:i+1000])+' '+\
            str(yearfile[:4])+'_'+('%06d'%(i,))[:3]+'.nc' 
            os.system(expr)

        #Once slices are done, combine
        expr = 'cdo -z zip -mergetime '+\
        str(yearfile[:4])+'_*.nc '+\
        str(yearfile[:4])+'.nc'
        os.system(expr)

        #Delete the slices
        expr = 'rm '+str(yearfile[:4])+'_*.nc'
        os.system(expr)


    else:
        expr = 'cdo -z zip -mergetime '+\
        str(yearfile[:4])+'/'+str(yearfile[:4])+'_*.nc '+\
        str(yearfile[:4])+'.nc' 
        os.system(expr)       

    #Move the remaining individual .nc files and the directory
    expr = 'rm -r pfiles_netcdf/'+str(yearfile[:4]) #Remove the old folder
    os.system(expr)
    expr = 'mv '+str(yearfile[:4]+' pfiles_netcdf')
    os.system(expr)

    expr = 'mv ' + yearfile + ' lists'
    os.system(expr)

    #Move the empty files list to the individual cast netcdf folder
    expr = 'mv emptyfile_problems.txt pfiles_netcdf/'+str(yearfile[:4])
    os.system(expr)
