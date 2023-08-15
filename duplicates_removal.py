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

Remove duplicates from the .pfile sourced data.
Duplicates are defined as two or more casts which occur within the same hour and are
less than 0.01 degrees distance from one another.

The cast which contains more temperature recordings is kept.

'''
#######################################################################################


#Create a list of files made with the outliers removed
path = '/gpfs/fs7/dfo/dpnm/joc000/Data/AZMP/Data_Input/NetCDF_Gen/empties_duplicates/'
nc_list = [path+'empties/'+i+'.nc' for i in np.arange(1912,2022+1).astype(str)]


#Cycle through each year
for i in np.sort(nc_list)[-28:]:

	#Load the .nc file, now with no empties
	pfile = xr.open_dataset(i)

	spatial_rounder = 2 #number of decimal points
	temporal_rounder = 'm' #minutes

	#Determine the times where there are time duplicates
	time = pfile.time.values.astype('datetime64['+temporal_rounder+']')
	longitude = pfile.longitude.values.round(spatial_rounder)
	latitude = pfile.latitude.values.round(spatial_rounder)

	#Round the time to the nearest 10th minute
	time_10min = np.zeros(time.size).astype(str)
	for ii,value in enumerate(time):
		minutes = int(str(value)[-1:])
		if minutes > 5:
			add = np.timedelta64(10-minutes,'m')
			time_10min[ii] = value+add
		else:
			subtract = np.timedelta64(minutes,'m')
			time_10min[ii] = value-subtract
	time_10min = np.array([np.datetime64(ii) for ii in time_10min])
	

	#Start to bring all the variables together
	spatiotemporal_comp = np.array([
		time_10min.astype(str),
		latitude.astype(str),
		longitude.astype(str)
		])

	a = spatiotemporal_comp.T

	#Merge the time,lat,lon as one string, consider as a whole
	a = [row[0]+','+row[1]+','+row[2] for row in spatiotemporal_comp.T]
	unique_dates,unique_index,unique_count = np.unique(a,return_counts=True,return_index=True)

	#Mark the times as duplicates
	duplicates_flag = np.full(time.shape,False) 
	
	for ii in np.arange(unique_dates.size):

		#Determine if there is a count higher than 1 (duplicate)
		if unique_count[ii] > 1:

			#Isolate for the casts
			spots = np.where(np.array(a) == unique_dates[ii])[0]


			#SPECIAL CASES
			#Determine if one of the casts is XBT and one is BO, if so skip
			if unique_count[ii] == 2:
				place1 = pfile.instrument_ID[spots].values
				if 'bo' in str(place1[0]).lower() or 'bo' in str(place1[1]).lower():
					if 'xb' in str(place1[0]).lower() or 'xb' in str(place1[1]).lower():
						continue

			#Determine if one of the casts is S1 and the other S0, and they're both fishsets
			if unique_count[ii] == 2:
				place1 = pfile.instrument_ID[spots].values
				place2 = pfile.comments[spots].values
				if 'fishset' in str(place2[0]).lower() and 'fishset' in str(place2[1]).lower():
					if 's0' in str(place1[0]).lower() or 's0' in str(place1[1]).lower():
						if 's1' in str(place1[0]).lower() or 's1' in str(place1[1]).lower():
	
							#Import the temperature
							temp = pfile.temperature[spots].values
							saln = pfile.salinity[spots].values

							#Determine the number of recordings in each cast
							nom_temp = (~np.isnan(temp)).sum(axis=1)
							nom_saln = (~np.isnan(saln)).sum(axis=1)							
							
							#Keep the cast with the highest number
							duplicates_flag[spots[np.argmin(nom_temp+nom_saln)]] = True
							continue

			#If repeat is the comment of each cast for two casts
			if unique_count[ii] == 2:
				place2 = pfile.comments[spots].values
				if 'repeat' in str(place2[0]).lower() and 'repeat' in str(place2[1]).lower():

					#Import the temperature
					temp = pfile.temperature[spots].values
					saln = pfile.salinity[spots].values

					#Determine the number of recordings in each cast
					nom_temp = (~np.isnan(temp)).sum(axis=1)
					nom_saln = (~np.isnan(saln)).sum(axis=1)							
					
					#Keep the cast with the highest number
					duplicates_flag[spots[np.argmin(nom_temp+nom_saln)]] = True
					continue

			#If repetitive is the comment of each cast for two casts
			if unique_count[ii] == 2:
				place2 = pfile.comments[spots].values
				if 'repetitive' in str(place2[0]).lower() and 'repetitive' in str(place2[1]).lower():

					#Import the temperature
					temp = pfile.temperature[spots].values
					saln = pfile.salinity[spots].values

					#Determine the number of recordings in each cast
					nom_temp = (~np.isnan(temp)).sum(axis=1)
					nom_saln = (~np.isnan(saln)).sum(axis=1)							
					
					#Keep the cast with the highest number
					duplicates_flag[spots[np.argmin(nom_temp+nom_saln)]] = True
					continue



			#Import the temperature
			temp = pfile.temperature[spots].values
			saln = pfile.salinity[spots].values

			#Determine the number of recordings in each cast
			nom_temp = (~np.isnan(temp)).sum(axis=1)
			nom_saln = (~np.isnan(saln)).sum(axis=1)

			#Plot the temperature, salinity
			plt.subplot(1,2,1)
			for iii in np.arange(temp.shape[0]):
				plt.scatter(temp[iii],np.arange(pfile.level.size),
					label='Temp'+str(iii)+': '+str(nom_temp[iii]))
			plt.legend()
			plt.gca().invert_yaxis()
			plt.xlim(xmin=-2,xmax=35)
			plt.grid()

			plt.subplot(1,2,2)
			for iii in np.arange(temp.shape[0]):
				plt.scatter(saln[iii],np.arange(pfile.level.size),
					label='Saln'+str(iii)+': '+str(nom_saln[iii]))
			plt.legend()
			plt.gca().invert_yaxis()
			plt.xlim(xmin=30,xmax=37)
			plt.grid()

			#Plot the cast info
			title = ''
			for iii in np.arange(temp.shape[0]):
				title = title+\
				'CAST '+str(iii)+': '+\
				pfile.instrument_ID[spots].values[iii]+', '+\
				pfile.comments[spots].values[iii]+', '+\
				pfile.file_names[spots].values[iii]+', '+\
				pfile.time[spots].values[iii].astype(str)+'\n'
			plt.suptitle(title,fontsize=8)

			#Determine the cast you would like to keep
			for test in range(1,10):
				user = input("Select profile number to keep (type 'all' to keep all).")
				while np.isin(user, np.concatenate((np.arange(spots.size),['all']))) == False:
					user = input('Not a valid choice. Please re-select a profile number.')
				if np.isin(user, np.arange(spots.size).astype(str)):
					#Record the kept cast
					x = int(user)
					duplicates_flag[spots[~np.isin(spots,spots[x])]] = True
					break
				elif user == 'all':
					break

			plt.close()
			print(str(ii)+' out of '+str(unique_dates.size)+' done.')


	#Remove the duplicates from the pfile
	pfile = pfile.sel(time=~duplicates_flag)

	#Flatten the values inside the string variables
	for ii in ['trip_ID','comments','instrument_ID','instrument_type','file_names']:
		pfile[ii] = pfile[ii].astype(str)


	#Save the new file
	pfile.to_netcdf(path+'duplicates/'+i[-7:-3]+'.nc',mode='w',\
		encoding={'time':{'units': "seconds since 1900-01-01 00:00:00",'calendar':	'gregorian'}})
	pfile.close()

	print(str(i[-7:-3])+' done.')
















