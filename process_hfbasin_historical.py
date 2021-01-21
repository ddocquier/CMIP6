#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Post-process northward ocean heat transport (hfbasin)
PROGRAMMER
    D. Docquier
LAST UPDATE
    28/08/2020
'''

# Standard libraries
from netCDF4 import Dataset
import numpy as np

# Options
institute = 'EC-Earth-Consortium'
model = 'EC-Earth3-Veg'
experiment = 'historical'
member = 'r6i1p1f1'
start_year = 1950
end_year = 2014
save_var = True # True: save variables in a .npy file; False: don't save variables

# Working directories
dir_input = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/' + str(member) + '/Omon/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/historical/hfbasin/'

# Load file for latitude dimension
if institute == 'EC-Earth-Consortium':
    file_init = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-195012.nc'
elif model == 'CNRM-CM6-1' or model == 'CNRM-ESM2-1':
    file_init = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/EC-Earth-Consortium/EC-Earth3/historical/r1i1p1f1/Omon/hfbasin/gn/latest/hfbasin_Omon_EC-Earth3_historical_r1i1p1f1_gn_195001-195012.nc'
elif model == 'MPI-ESM1-2-LR' or model == 'MPI-ESM1-2-HAM':
    file_init = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-196912.nc'
elif model == 'MPI-ESM1-2-HR':
    file_init = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-195412.nc'
elif institute == 'MRI':
    file_init = dir_input + 'hfbasin/gr2z/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr2z_185001-201412.nc'
elif model == 'NorESM2-LM' or model == 'NorESM2-MM':
    file_init = dir_input + 'hfbasin/grz/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_grz_195001-195912.nc'
elif model == 'NorCPM1':
    file_init = dir_input + 'hfbasin/grz/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_grz_185001-201412.nc'
elif model == 'SAM0-UNICON':
    file_init = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-195912.nc'
elif model == 'HadGEM3-GC31-LL' or model == 'UKESM1-0-LL':
    file_init = dir_input + 'hfbasin/gnz/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gnz_195001-201412.nc'
elif model == 'HadGEM3-GC31-MM':
    file_init = dir_input + 'hfbasin/gnz/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gnz_195001-196912.nc'
elif institute == 'NASA-GISS':
    file_init = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195101-200012.nc'
else:
    file_init = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_185001-201412.nc'    
fh = Dataset(file_init, mode='r')
if institute == 'IPSL':
    lat = fh.variables['nav_lat'][:,0]
    nlat = np.size(lat)
else:
    lat = fh.variables['lat'][:]
    nlat = np.size(lat)
fh.close()

# Initialization
nyears = int(end_year-start_year+1)
nmy = int(12)
nm = int(nyears * nmy)
hfbasin_total = np.zeros((nyears,nlat))
hfbasin_atlantic = np.zeros((nyears,nlat))
hfbasin_pacific = np.zeros((nyears,nlat))

# Load and store variables
if institute == 'EC-Earth-Consortium':
    
    # Loop over years
    for year in np.arange(nyears):
        print(start_year+year)

        # Retrieve hfbasin
        if (member == 'r10i1p1f1' and (start_year+year) == 1996) or (member == 'r21i1p1f1' and (start_year+year) == 1969):
            file_hfbasin = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/hfbasin/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
        else:
            file_hfbasin = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
        fh = Dataset(file_hfbasin, mode='r')
        hfbasin_tot = fh.variables['hfbasin'][:,0,:]
        hfbasin_atl = fh.variables['hfbasin'][:,1,:]
        hfbasin_pac = fh.variables['hfbasin'][:,2,:]
        fh.close()

        # Compute annual mean and store variables
        hfbasin_total[year,:] = np.nanmean(hfbasin_tot/1.e15,axis=0)
        hfbasin_atlantic[year,:] = np.nanmean(hfbasin_atl/1.e15,axis=0)
        hfbasin_pacific[year,:] = np.nanmean(hfbasin_pac/1.e15,axis=0)

elif institute == 'MPI-M' or institute == 'HAMMOZ-Consortium' or model == 'NorESM2-LM' or model == 'NorESM2-MM' or model == 'SAM0-UNICON' or model == 'HadGEM3-GC31-MM' or institute == 'NASA-GISS':
    
    # Loop over files
    if model == 'MPI-ESM1-2-HR':
        n_files = 13
    elif model == 'NorESM2-LM' or model == 'NorESM2-MM' or model == 'SAM0-UNICON':
        n_files = 7
    elif model == 'HadGEM3-GC31-MM':
        n_files = 5
    elif institute == 'NASA-GISS':
        n_files = 3
    else:
        n_files = 4

    for filex in np.arange(n_files):
        print('file',filex+1)
        if model == 'MPI-ESM1-2-HR':
            nyears_file = 5
            j = nyears_file - 1
        elif model == 'NorESM2-LM' or model == 'NorESM2-MM' or model == 'SAM0-UNICON':
            nyears_file = 10
            if filex == 6:
                j = 4
            else:
                j = nyears_file - 1
        elif model == 'HadGEM3-GC31-MM':
            if filex == 0:
                start_period = '1950'
                end_period = '1969'
                nyears_file = 20
            elif filex == 1:
                start_period = '1970'
                end_period = '1987'
                nyears_file = 18
            elif filex == 2:
                start_period = '1988'
                end_period = '1989'
                nyears_file = 2
            elif filex == 3:
                start_period = '1990'
                end_period = '2009'
                nyears_file = 20
            elif filex == 4:
                start_period = '2010'
                end_period = '2014'
                nyears_file = 5
        elif institute == 'NASA-GISS':
            if filex == 0:
                start_period = '1901'
                end_period = '1950'
                nyears_file = 1
            elif filex == 1:
                start_period = '1951'
                end_period = '2000'
                nyears_file = 50
            elif filex == 2:
                start_period = '2001'
                end_period = '2014'
        else:
            nyears_file = 20
            if filex == 3:
                j = 4
            else:
                j = nyears_file - 1
    
        # Retrieve hfbasin
        if model == 'NorESM2-LM' or model == 'NorESM2-MM':
            file_hfbasin = dir_input + 'hfbasin/grz/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_grz_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
        elif model == 'HadGEM3-GC31-MM':
            file_hfbasin = dir_input + 'hfbasin/gnz/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gnz_' + str(start_period) + '01-' + str(end_period) + '12.nc'
        elif institute == 'NASA-GISS':
            file_hfbasin = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_period) + '01-' + str(end_period) + '12.nc'
        else:
            file_hfbasin = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
        fh = Dataset(file_hfbasin, mode='r')
        if model == 'NorESM2-LM' or model == 'NorESM2-MM':
            hfbasin_atl = fh.variables['hfbasin'][:,1,:]
            hfbasin_pac = fh.variables['hfbasin'][:,2,:]
            hfbasin_tot = fh.variables['hfbasin'][:,3,:]
        elif model == 'SAM0-UNICON' or institute == 'NASA-GISS':
            hfbasin_atl = fh.variables['hfbasin'][:,0,:]
            hfbasin_pac = fh.variables['hfbasin'][:,1,:]
            hfbasin_tot = fh.variables['hfbasin'][:,2,:]
        elif model == 'HadGEM3-GC31-MM':
            hfbasin_atl = fh.variables['hfbasin'][:,0,:]
            hfbasin_tot = fh.variables['hfbasin'][:,1,:]
            hfbasin_pac = fh.variables['hfbasin'][:,2,:]
        else:
            hfbasin_tot = fh.variables['hfbasin'][:,0,:]
            hfbasin_atl = fh.variables['hfbasin'][:,1,:]
            hfbasin_pac = fh.variables['hfbasin'][:,2,:]
        nt,notused = hfbasin_atl.shape
        fh.close()

        # Store variables
        subtot_yrs = int(nt / nmy)
        for year in np.arange(subtot_yrs):
            if model == 'HadGEM3-GC31-MM':
                if filex == 0:
                    hfbasin_total[year,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_atlantic[year,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_pacific[year,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                elif filex == 1:
                    hfbasin_total[year+20,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_atlantic[year+20,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_pacific[year+20,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                elif filex == 2:
                    hfbasin_total[year+38,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_atlantic[year+38,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_pacific[year+38,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                elif filex == 3:
                    hfbasin_total[year+40,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_atlantic[year+40,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_pacific[year+40,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                elif filex == 4:
                    hfbasin_total[year+60,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_atlantic[year+60,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_pacific[year+60,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
            elif institute == 'NASA-GISS':
                if filex == 0:
                    hfbasin_total[year,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_atlantic[year,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_pacific[year,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                elif filex == 1:
                    hfbasin_total[year+1,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_atlantic[year+1,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_pacific[year+1,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                elif filex == 2:
                    hfbasin_total[year+51,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_atlantic[year+51,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                    hfbasin_pacific[year+51,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
            else:
                hfbasin_total[year+filex*nyears_file,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                hfbasin_atlantic[year+filex*nyears_file,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
                hfbasin_pacific[year+filex*nyears_file,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)

else:

    # Retrieve hfbasin
    if institute == 'MRI':
        file_hfbasin = dir_input + 'hfbasin/gr2z/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr2z_185001-201412.nc'
    elif model == 'NorCPM1':
        file_hfbasin = dir_input + 'hfbasin/grz/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_grz_185001-201412.nc'
    elif model == 'HadGEM3-GC31-LL' or model == 'UKESM1-0-LL':
        file_hfbasin = dir_input + 'hfbasin/gnz/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gnz_195001-201412.nc'
    else:
        file_hfbasin = dir_input + 'hfbasin/gn/latest/hfbasin_Omon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_185001-201412.nc'
    fh = Dataset(file_hfbasin, mode='r')
    if institute == 'CAS':
        hfbasin_tot = fh.variables['hfbasin'][:,0,:]
        hfbasin_atl = fh.variables['hfbasin'][:,1,:]
        hfbasin_pac = hfbasin_tot - hfbasin_atl
    elif model == 'CNRM-CM6-1' or model == 'CNRM-ESM2-1':
        hfbasin_tot = fh.variables['hfbasin'][:,0,1:-1]
        hfbasin_atl = fh.variables['hfbasin'][:,1,1:-1]
        hfbasin_pac = fh.variables['hfbasin'][:,2,1:-1]
    elif institute == 'CCCma' or institute == 'MRI' or model == 'NorCPM1':
        hfbasin_atl = fh.variables['hfbasin'][:,0,:]
        hfbasin_pac = fh.variables['hfbasin'][:,1,:]
        hfbasin_tot = fh.variables['hfbasin'][:,2,:]
    elif institute == 'MOHC':
        hfbasin_atl = fh.variables['hfbasin'][:,0,:]
        hfbasin_tot = fh.variables['hfbasin'][:,1,:]
        hfbasin_pac = fh.variables['hfbasin'][:,2,:]
    elif institute == 'IPSL':
        hfbasin_tot = fh.variables['hfbasin'][:,0,:,0]
        hfbasin_atl = fh.variables['hfbasin'][:,1,:,0]
        hfbasin_pac = fh.variables['hfbasin'][:,2,:,0]
    fh.close()

    # Compute annual mean and store variables
    for year in np.arange(nyears):
        hfbasin_total[year,:] = np.nanmean(hfbasin_tot[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
        hfbasin_atlantic[year,:] = np.nanmean(hfbasin_atl[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)
        hfbasin_pacific[year,:] = np.nanmean(hfbasin_pac[year*nmy:year*nmy+nmy-1,:]/1.e15,axis=0)


# Save variables
if save_var == True:
    filename = dir_output + 'hfbasin_annmean_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year) + '-' +str(end_year) + '.npy'
    np.save(filename,[hfbasin_total,hfbasin_atlantic,hfbasin_pacific,lat])
