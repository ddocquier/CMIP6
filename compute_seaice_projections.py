#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Compute monthly mean Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3)
    CMIP6 projections
PROGRAMMER
    D. Docquier
LAST UPDAT9
    06/10/2020
'''

# Standard libraries
from netCDF4 import Dataset
import numpy as np

# Options
institute = 'MIROC'
model = 'MIROC6'
experiment = 'ssp126' # ssp126; ssp585
member = 'r50i1p1f1'
start_year = 2015 # 2015
end_year = 2100 # CAMS: 2099; else: 2100
save_var = True # True: save variables in a .npy file; False: don't save variables

# Parameters
lat_threshold = 40.

# Working directories
if institute == 'AWI' or model == 'KIOST-ESM' or model == 'MIROC6':
    dir_input = '/home/rossby/data_lib/esgf/CMIP6/ScenarioMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/' + str(member) + '/SImon/'
else:
    dir_input = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/ScenarioMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/' + str(member) + '/SImon/'
if institute == 'CNRM-CERFACS' or model == 'MIROC-ES2L':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/historical/r1i1p1f2/Ofx/areacello/gn/latest/'
elif model == 'CanESM5-CanOE':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/historical/r1i1p2f1/Ofx/areacello/gn/latest/'
elif model == 'CESM2-WACCM':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/CESM2/historical/r1i1p1f1/Ofx/areacello/gn/latest/'
elif model == 'FGOALS-g3':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/FGOALS-f3-L/historical/r1i1p1f1/Ofx/areacello/gn/latest/'
elif model == 'TaiESM1' or model == 'CAMS-CSM1-0' or institute == 'E3SM-Project' or model == 'FIO-ESM-2-0' or institute == 'INM' or institute == 'MOHC' or model == 'GISS-E2-1-H' or model == 'NESM3' or institute == 'BCC' or model == 'KIOST-ESM':
    dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/areacello/'
    if model == 'GISS-E2-1-H':
        dir_grid2 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/areacella/'
elif model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
    dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/areacella/'
elif model == 'MPI-ESM1-2-HR':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/MPI-M/' + str(model) + '/historical/r1i1p1f1/Ofx/areacello/gn/latest/'
else:
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/historical/r1i1p1f1/Ofx/areacello/gn/latest/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/' + str(experiment) + '/SIA-SIV/'

# Function to compute total Arctic sea-ice area (10^6 km^2)
def compute_area(nm,lat,lat_threshold,siconc,grid_area):
    area = np.zeros(nm)
    for i in np.arange(nm):
        area[i] = np.nansum((lat >= lat_threshold) * (siconc[i,:,:] / 100.) * grid_area)
    area = area / 1.e6
    area[area <= 0.] = np.nan
    return area

# Function to compute total Arctic sea-ice volume (10^3 km^3)
def compute_vol(nm,lat,lat_threshold,sivol,grid_area):
    volume = np.zeros(nm)
    for i in np.arange(nm):
        volume[i] = np.nansum((lat >= lat_threshold) * (sivol[i,:,:] / 1.e3) * grid_area)
    volume = volume / 1.e3
    volume[volume <= 0.] = np.nan
    return volume

# Load grid-cell area (km^2) and latitude
if not (institute == 'AWI'):
    if institute == 'CNRM-CERFACS' or model == 'MIROC-ES2L':
        file_grid = dir_grid + 'areacello_Ofx_' + str(model) + '_historical_r1i1p1f2_gn.nc'
    elif model == 'CanESM5-CanOE':
        file_grid = dir_grid + 'areacello_Ofx_' + str(model) + '_historical_r1i1p2f1_gn.nc'
    elif model == 'CESM2-WACCM':
        file_grid = dir_grid + 'areacello_Ofx_CESM2_historical_r1i1p1f1_gn.nc'
    elif model == 'TaiESM1' or model == 'CAMS-CSM1-0' or institute == 'E3SM-Project' or model == 'FIO-ESM-2-0' or institute == 'INM' or institute == 'MOHC' or model == 'GISS-E2-1-H' or model == 'NESM3' or model == 'KIOST-ESM':
        file_grid = dir_grid + 'areacello_' + str(model) + '.nc'
        if model == 'GISS-E2-1-H':
            file_grid2 = dir_grid2 + 'areacella_GISS-E2-1-G.nc'
    elif model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
        file_grid = dir_grid + 'areacella_GISS-E2-1-G.nc'
    elif institute == 'BCC':
        file_grid = dir_grid + 'areacello_' + str(institute) + '.nc'
    elif model == 'FGOALS-g3':
        file_grid = dir_grid + 'areacello_Ofx_FGOALS-f3-L_historical_r1i1p1f1_gn.nc'
    else:
        file_grid = dir_grid + 'areacello_Ofx_' + str(model) + '_historical_r1i1p1f1_gn.nc'
    fh = Dataset(file_grid, mode='r')
    if model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
        grid_area = fh.variables['areacella'][:]
    else:
        grid_area = fh.variables['areacello'][:]
    if model == 'GISS-E2-1-H':
        fh2 = Dataset(file_grid2, mode='r')
        grid_area2 = fh2.variables['areacella'][:]
        grid_area2 = grid_area2 / 1.e6
    if model == 'NorESM2-LM' or model == 'NorESM2-MM':
        grid_area = grid_area[0:-1,:]
    grid_area = grid_area / 1.e6
    if institute == 'CNRM-CERFACS' or institute == 'NCAR' or institute == 'NOAA-GFDL' or model == 'NESM3':
        lat = fh.variables['lat'][:]
    elif institute == 'IPSL':
        lat = fh.variables['nav_lat'][:]
    elif institute == 'E3SM-Project' or institute == 'INM' or institute == 'NASA-GISS' or model == 'KIOST-ESM':
        lat_init = fh.variables['lat'][:]
        lon_init = fh.variables['lon'][:]
        lon,lat = np.meshgrid(lon_init,lat_init)
        if model == 'GISS-E2-1-H':
            lat_init2 = fh2.variables['lat'][:]
            lon_init2 = fh2.variables['lon'][:]
            lon2,lat2 = np.meshgrid(lon_init2,lat_init2)
    else:
        lat = fh.variables['latitude'][:]
        if institute == 'CAS':
            grid_area = np.flipud(grid_area)
            lat = np.flipud(lat)
        if model == 'NorESM2-LM' or model == 'NorESM2-MM':
            lat = lat[0:-1,:]
    fh.close()
    if model == 'GISS-E2-1-H':
        fh2.close()

# Initialization
nyears = int(end_year-start_year+1)
nmy = int(12)
nm = int(nyears * nmy)
area_total = np.zeros((nyears,nmy))
volume_total = np.zeros((nyears,nmy))

# Load and store variables
if institute == 'EC-Earth-Consortium':
    
    # Loop over years
    for year in np.arange(nyears):
        print(start_year+year)

        # Retrieve siconc
        if model == 'EC-Earth3' and experiment == 'ssp585' and member == 'r9i1p1f1':
            file_siconc = dir_input + 'siconc/gn/v20190710/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
        else:
            file_siconc = dir_input + 'siconc/gn/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
        fh = Dataset(file_siconc, mode='r')
        siconc = fh.variables['siconc'][:]
        siconc[siconc < 0.] = np.nan
        siconc[siconc > 101.] = np.nan
        fh.close()

        # Retrieve sivol
        #file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
        #fh = Dataset(file_sivol, mode='r')
        #sivol = fh.variables['sivol'][:]
        #sivol[sivol < 0.] = np.nan
        #sivol[sivol > 20.] = np.nan
        #fh.close()

        # Retrieve sithick and compute sivol
        if model == 'EC-Earth3' and experiment == 'ssp585' and member == 'r9i1p1f1':
            sivol = siconc * np.nan
        else:
            file_sithick = '/home/rossby/data_lib/esgf/CMIP6/ScenarioMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/' + str(member) + '/SImon/sithick/gn/latest/sithick_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
            fh = Dataset(file_sithick, mode='r')
            sithick = fh.variables['sithick'][:]
            sithick[sithick < 0.] = np.nan
            sithick[sithick > 40.] = np.nan
            fh.close()
            sivol = (sithick * siconc) / 100.

        # Compute total Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3)
        area_month = compute_area(nmy,lat,lat_threshold,siconc,grid_area)
        volume_month = compute_vol(nmy,lat,lat_threshold,sivol,grid_area)

        # Store variables
        for month in np.arange(nmy):
            area_total[year,month] = area_month[month]
            volume_total[year,month] = volume_month[month]

elif model == 'MPI-ESM1-2-HR' or institute == 'MPI-M' or model == 'CESM2' or (model == 'CESM2-WACCM' and member == 'r3i1p1f1') or (model == 'CESM2-WACCM' and member == 'r4i1p1f1') or (model == 'CESM2-WACCM' and member == 'r5i1p1f1') or model == 'NorESM2-LM' or model == 'NorESM2-MM' or model == 'SAM0-UNICON' or institute == 'E3SM-Project' or institute == 'INM' or institute == 'MOHC' or institute == 'NASA-GISS' or institute == 'AWI':

    # Loop over files
    if model == 'MPI-ESM1-2-HR' or institute == 'E3SM-Project':
        n_files = 18
    elif model == 'NorESM2-LM' or model == 'NorESM2-MM' or institute == 'AWI':
        n_files = 9
    elif institute == 'NCAR' or institute == 'INM' or model == 'HadGEM3-GC31-LL' or model == 'UKESM1-0-LL':
        n_files = 2
    elif model == 'HadGEM3-GC31-MM':
        n_files = 5
    elif model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
        n_files = 3
    else:
        n_files = 5

    for filex in np.arange(n_files):
        print('file',filex+1)
        if model == 'MPI-ESM1-2-HR' or institute == 'E3SM-Project':
            nyears_file = 5
            if filex == 17:
                j = 0
            else:
                j = nyears_file - 1
        elif model == 'NorESM2-LM' or model == 'NorESM2-MM' or institute == 'AWI':
            if filex == 0:
                start_period = '2015'
                end_period = '2020'
                nyears_file = 6
            elif filex == 1:
                start_period = '2021'
                end_period = '2030'
                nyears_file = 10
            elif filex == 2:
                start_period = '2031'
                end_period = '2040'
                nyears_file = 10
            elif filex == 3:
                start_period = '2041'
                end_period = '2050'
                nyears_file = 10
            elif filex == 4:
                start_period = '2051'
                end_period = '2060'
                nyears_file = 10
            elif filex == 5:
                start_period = '2061'
                end_period = '2070'
                nyears_file = 10
            elif filex == 6:
                start_period = '2071'
                end_period = '2080'
                nyears_file = 10
            elif filex == 7:
                start_period = '2081'
                end_period = '2090'
                nyears_file = 10
            elif filex == 8:
                start_period = '2091'
                end_period = '2100'
                nyears_file = 10
        elif institute == 'NCAR' or institute == 'INM':
            nyears_file = 50
            if filex == 0:
                j = nyears_file - 1
            elif filex == 1:
                j = 35
        elif model == 'HadGEM3-GC31-LL' or model == 'UKESM1-0-LL':
            nyears_file = 35
            if filex == 0:
                j = nyears_file - 1
            elif filex == 1:
                j = 50
        elif model == 'HadGEM3-GC31-MM':
            if filex == 0:
                start_period = '2015'
                end_period = '2029'
                nyears_file = 15
            elif filex == 1:
                start_period = '2030'
                end_period = '2049'
                nyears_file = 20
            elif filex == 2:
                start_period = '2050'
                end_period = '2069'
                nyears_file = 20
            elif filex == 3:
                start_period = '2070'
                end_period = '2089'
                nyears_file = 20
            elif filex == 4:
                start_period = '2090'
                end_period = '2100'
                nyears_file = 11
        elif model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
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
            if filex == 4:
                j = 5
            else:
                j = nyears_file - 1

        # Retrieve siconc
        if not (institute == 'AWI' or model == 'GISS-E2-1-G-CC'):
            if institute == 'E3SM-Project' or model == 'GISS-E2-1-H':
                file_siconc = dir_input + 'siconc/gr/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
            elif institute == 'INM':
                file_siconc = dir_input + 'siconc/gr1/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr1_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
            elif model == 'HadGEM3-GC31-MM' or institute == 'NCC':
                file_siconc = dir_input + 'siconc/gn/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_period) + '01-' + str(end_period) + '12.nc'
            elif model == 'GISS-E2-1-G':
                file_siconc = dir_input + 'siconca/gn/latest/siconca_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_period) + '01-' + str(end_period) + '12.nc'
            else:
                file_siconc = dir_input + 'siconc/gn/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
            fh = Dataset(file_siconc, mode='r')
            if model == 'GISS-E2-1-G':
                siconc = fh.variables['siconca'][:]
            else:
                siconc = fh.variables['siconc'][:]
            if model == 'GISS-E2-1-G' and filex == 0:
                siconc = siconc[588::,:,:]
            siconc[siconc < 0.] = np.nan
            siconc[siconc > 101.] = np.nan
            nt,notused,notused = siconc.shape
            fh.close()

        # Retrieve sivol
        if not (institute == 'E3SM-Project' or institute == 'INM' or institute == 'AWI' or model == 'GISS-E2-1-H' or (model == 'CESM2' and member == 'r1i1p1f1') or (model == 'CESM2' and member == 'r2i1p1f1')):
            if model == 'HadGEM3-GC31-MM' or model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC' or institute == 'NCC':
                file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_period) + '01-' + str(end_period) + '12.nc'
            elif institute == 'AWI':
                file_sivol = dir_input + 'sivoln/gn/latest/sivoln_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_period) + '01-' + str(end_period) + '12.nc'
            else:
                file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
            fh = Dataset(file_sivol, mode='r')
            sivol = fh.variables['sivol'][:]
            if institute == 'AWI':
                sivol = sivol[:,0]
                if filex == 0:
                    sivol = sivol[108::]
            if (model == 'GISS-E2-1-G' and filex == 0) or (model == 'GISS-E2-1-G-CC' and filex == 0):
                sivol = sivol[588::,:,:]
            sivol[sivol < 0.] = np.nan
            sivol[sivol > 20.] = np.nan
            if model == 'GISS-E2-1-G-CC':
                nt,notused,notused = sivol.shape
            fh.close()

        # Retrieve siarean
        if institute == 'AWI':
            file_siarean = dir_input + 'siarean/gn/latest/siarean_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_period) + '01-' + str(end_period) + '12.nc'
            fh = Dataset(file_siarean, mode='r')
            area_month = fh.variables['siarean'][:]
            area_month = area_month[:,0]
            area_month[area_month < 0.] = np.nan
            area_month[area_month > 40.] = np.nan
            nt = np.size(area_month)
            fh.close()

        # Retrieve sivoln
        if institute == 'AWI':
            file_sivoln = dir_input + 'sivoln/gn/latest/sivoln_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_period) + '01-' + str(end_period) + '12.nc'
            fh = Dataset(file_sivoln, mode='r')
            volume_month = fh.variables['sivoln'][:]
            volume_month = volume_month[:,0]
            volume_month = volume_month / 1.e3
            volume_month[volume_month < 0.] = np.nan
            volume_month[volume_month > 100.] = np.nan
            fh.close()

        # Retrieve sithick and compute sivol
        if institute == 'E3SM-Project':
            file_sithick = dir_input + 'sithick/gr/latest/sithick_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
            fh = Dataset(file_sithick, mode='r')
            sithick = fh.variables['sithick'][:]
            sithick[sithick < 0.] = np.nan
            sithick[sithick > 40.] = np.nan
            fh.close()
            sivol = (sithick * siconc) / 100.

        # Compute total Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3)
        if not (institute == 'AWI' or model == 'GISS-E2-1-G-CC'):
            area_month = compute_area(nt,lat,lat_threshold,siconc,grid_area)
        if not (institute == 'INM' or model == 'GISS-E2-1-H' or institute == 'AWI' or (model == 'CESM2' and member == 'r1i1p1f1') or (model == 'CESM2' and member == 'r2i1p1f1')):
            volume_month = compute_vol(nt,lat,lat_threshold,sivol,grid_area)

        # Store variables
        subtot_yrs = int(nt / nmy)
        for year in np.arange(subtot_yrs):
            for month in np.arange(nmy):
                if model == 'HadGEM3-GC31-MM':
                    if filex == 0:
                        area_total[year,month] = area_month[month+year*nmy]
                        volume_total[year,month] = volume_month[month+year*nmy]
                    elif filex == 1:
                        area_total[year+15,month] = area_month[month+year*nmy]
                        volume_total[year+15,month] = volume_month[month+year*nmy]
                    elif filex == 2:
                        area_total[year+35,month] = area_month[month+year*nmy]
                        volume_total[year+35,month] = volume_month[month+year*nmy]
                    elif filex == 3:
                        area_total[year+55,month] = area_month[month+year*nmy]
                        volume_total[year+55,month] = volume_month[month+year*nmy]
                    elif filex == 4:
                        area_total[year+75,month] = area_month[month+year*nmy]
                        volume_total[year+75,month] = volume_month[month+year*nmy]
                elif institute == 'NCC' or institute == 'AWI':
                    if filex == 0:
                        area_total[year,month] = area_month[month+year*nmy]
                        volume_total[year,month] = volume_month[month+year*nmy]
                    elif filex == 1:
                        area_total[year+6,month] = area_month[month+year*nmy]
                        volume_total[year+6,month] = volume_month[month+year*nmy]
                    elif filex == 2:
                        area_total[year+16,month] = area_month[month+year*nmy]
                        volume_total[year+16,month] = volume_month[month+year*nmy]
                    elif filex == 3:
                        area_total[year+26,month] = area_month[month+year*nmy]
                        volume_total[year+26,month] = volume_month[month+year*nmy]
                    elif filex == 4:
                        area_total[year+36,month] = area_month[month+year*nmy]
                        volume_total[year+36,month] = volume_month[month+year*nmy]
                    elif filex == 5:
                        area_total[year+46,month] = area_month[month+year*nmy]
                        volume_total[year+46,month] = volume_month[month+year*nmy]
                    elif filex == 6:
                        area_total[year+56,month] = area_month[month+year*nmy]
                        volume_total[year+56,month] = volume_month[month+year*nmy]
                    elif filex == 7:
                        area_total[year+66,month] = area_month[month+year*nmy]
                        volume_total[year+66,month] = volume_month[month+year*nmy]
                    elif filex == 8:
                        area_total[year+76,month] = area_month[month+year*nmy]
                        volume_total[year+76,month] = volume_month[month+year*nmy]
                elif model == 'GISS-E2-1-G':
                    if filex == 0:
                        area_total[year,month] = area_month[month+year*nmy]
                        volume_total[year,month] = volume_month[month+year*nmy]
                    elif filex == 1:
                        area_total[year+1,month] = area_month[month+year*nmy]
                        volume_total[year+1,month] = volume_month[month+year*nmy]
                    elif filex == 2:
                        area_total[year+51,month] = area_month[month+year*nmy]
                        volume_total[year+51,month] = volume_month[month+year*nmy]
                elif model == 'GISS-E2-1-G-CC':
                    if filex == 0:
                        area_total[year,month] = np.nan
                        volume_total[year,month] = volume_month[month+year*nmy]
                    elif filex == 1:
                        area_total[year+1,month] = np.nan
                        volume_total[year+1,month] = volume_month[month+year*nmy]
                    elif filex == 2:
                        area_total[year+51,month] = np.nan
                        volume_total[year+51,month] = volume_month[month+year*nmy]
                else:
                    area_total[year+filex*nyears_file,month] = area_month[month+year*nmy]
                    if institute == 'INM' or model == 'GISS-E2-1-H' or (model == 'CESM2' and member == 'r1i1p1f1') or (model == 'CESM2' and member == 'r2i1p1f1'):
                        volume_total[year,month] = np.nan
                    else:
                        volume_total[year+filex*nyears_file,month] = volume_month[month+year*nmy]

    # Loop over sivol files
    if model == 'GISS-E2-1-H':
        n_files = 3
        for filex in np.arange(n_files):
            print('file',filex+1)
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
            
            # Retrieve sivol
            file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_period) + '01-' + str(end_period) + '12.nc'
            fh = Dataset(file_sivol, mode='r')
            sivol = fh.variables['sivol'][:]
            if filex == 0:
                sivol = sivol[588::,:,:]
            sivol[sivol < 0.] = np.nan
            sivol[sivol > 20.] = np.nan
            nt,notused,notused = sivol.shape
            fh.close()

            # Compute total Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3)
            volume_month = compute_vol(nt,lat2,lat_threshold,sivol,grid_area2)

            # Store variables
            subtot_yrs = int(nt / nmy)
            for year in np.arange(subtot_yrs):
                for month in np.arange(nmy):
                    if filex == 0:
                        volume_total[year,month] = volume_month[month+year*nmy]
                    elif filex == 1:
                        volume_total[year+1,month] = volume_month[month+year*nmy]
                    elif filex == 2:
                        volume_total[year+51,month] = volume_month[month+year*nmy]

else:

    # Retrieve siconc
    if institute == 'CAMS':
        file_siconc = dir_input + 'siconc/gn/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_201501-209912.nc'
    else:
        if (model == 'CanESM5' and experiment == 'ssp585' and member == 'r17i1p2f1'):
            file_siconc = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/ssp585/siconc/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_201501-210012.nc'
        elif model == 'FIO-ESM-2-0':
            file_siconc = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/' + str(experiment) + '/siconc/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_201501-210012.nc'
        elif model == 'KIOST-ESM':
            file_siconc = dir_input + 'siconc/gr1/v20200911/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr1_201501-210012.nc'
        else:
            file_siconc = dir_input + 'siconc/gn/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_201501-210012.nc'
    fh = Dataset(file_siconc, mode='r')
    siconc = fh.variables['siconc'][:]
    siconc[siconc < 0.] = np.nan
    siconc[siconc > 101.] = np.nan
    fh.close()

    # Retrieve sivol
    if not (institute == 'CCCma' or institute == 'MIROC' or model == 'NorCPM1' or model == 'FGOALS-g3' or model == 'NESM3' or (model == 'GFDL-ESM4' and experiment == 'ssp585') or model == 'KIOST-ESM'):
        if institute == 'CAMS':
            file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_201501-209912.nc'
        elif model == 'FIO-ESM-2-0':
            file_sivol = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/' + str(experiment) + '/sivol/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_201501-210012.nc'
        else:
            file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_201501-210012.nc'
        fh = Dataset(file_sivol, mode='r')
        sivol = fh.variables['sivol'][:]
        sivol[sivol < 0.] = np.nan
        sivol[sivol > 20.] = np.nan
        fh.close()

    # Retrieve sithick and compute sivol
    if model == 'CanESM5' or institute == 'MIROC' or model == 'NESM3' or model == 'NorCPM1' or model == 'KIOST-ESM':
       if model == 'KIOST-ESM':
           file_sithick = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/' + str(experiment) + '/sithick/sithick_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr1_201501-210012.nc'
       else:
           file_sithick = '/home/rossby/data_lib/esgf/CMIP6/ScenarioMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/' + str(member) + '/SImon/sithick/gn/latest/sithick_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_201501-210012.nc'
       fh = Dataset(file_sithick, mode='r')
       sithick = fh.variables['sithick'][:]
       if (model == 'KIOST-ESM' and experiment == 'ssp126'):
           sithick = np.insert(sithick,108,np.nan,axis=0) # 01/2024
           sithick = np.insert(sithick,109,np.nan,axis=0) # 02/2024
           sithick = np.insert(sithick,110,np.nan,axis=0) # 03/2024
           sithick = np.insert(sithick,111,np.nan,axis=0) # 04/2024
           sithick = np.insert(sithick,112,np.nan,axis=0) # 05/2024
           sithick = np.insert(sithick,113,np.nan,axis=0) # 06/2024
           sithick = np.insert(sithick,114,np.nan,axis=0) # 07/2024
           sithick = np.insert(sithick,115,np.nan,axis=0) # 08/2024
           sithick = np.insert(sithick,116,np.nan,axis=0) # 09/2024
           sithick = np.insert(sithick,117,np.nan,axis=0) # 10/2024
           sithick = np.insert(sithick,118,np.nan,axis=0) # 11/2024
           sithick = np.insert(sithick,119,np.nan,axis=0) # 12/2024
       sithick[sithick < 0.] = np.nan
       sithick[sithick > 40.] = np.nan
       fh.close()
       sivol = (sithick * siconc) / 100.

    # Compute total Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3)
    area_month = compute_area(nm,lat,lat_threshold,siconc,grid_area)
    if not (model == 'FGOALS-g3' or (model == 'GFDL-ESM4' and experiment == 'ssp585')):
        volume_month = compute_vol(nm,lat,lat_threshold,sivol,grid_area)

    # Store variables
    for month in np.arange(nmy):
        for year in np.arange(nyears):
            area_total[year,month] = area_month[month+year*nmy]
            if model == 'FGOALS-g3' or (model == 'GFDL-ESM4' and experiment == 'ssp585'):
                volume_total[year,month] = np.nan
            else:
                volume_total[year,month] = volume_month[month+year*nmy]

    # Retrieve sea-ice area and volume for models that do not provide siconc and sivol
    if model == 'CIESM':
        file_siarean = dir_input + 'siarean/gn/latest/siarean_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-201412.nc'
        fh = Dataset(file_siarean, mode='r')
        siarean = fh.variables['siarean'][:]
        siarean[siarean < 0.] = np.nan
        siarean[siarean > 40.] = np.nan
        fh.close()
        file_sivoln = dir_input + 'sivoln/gn/latest/sivoln_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-201412.nc'
        fh = Dataset(file_sivoln, mode='r')
        sivoln = fh.variables['sivoln'][:]
        sivoln[sivoln < 0.] = np.nan
        sivoln[sivoln > 100.] = np.nan
        fh.close()
        print(siarean[8:nm:12])
        print(sivoln[8:nm:12])
        area_total = np.zeros((nyears,nmy))
        volume_total = np.zeros((nyears,nmy))
        for i in np.arange(nmy):
            for k in np.arange(nyears):
                area_total[k,i] = siarean[i+k*nmy]
                volume_total[k,i] = sivoln[i+k*nmy]

# Save variables
if save_var == True:
    filename = dir_output + 'SIA-SIV_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year) + '-' +str(end_year) + '.npy'
    np.save(filename,[area_total,volume_total])
