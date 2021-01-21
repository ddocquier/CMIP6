#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Compute monthly mean Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3)
    CMIP6 historical
PROGRAMMER
    D. Docquier
LAST UPDATE
    06/10/2020
'''

# Standard libraries
from netCDF4 import Dataset
import numpy as np

# Options
institute = 'KIOST'
model = 'KIOST-ESM'
experiment = 'historical'
member = 'r1i1p1f1'
start_year = 1950
end_year = 2014
save_var = True # True: save variables in a .npy file; False: don't save variables

# Parameters
lat_threshold = 40.

# Working directories
if model == 'CMCC-CM2-HR4' or model == 'KIOST-ESM':
    dir_input = '/home/rossby/data_lib/esgf/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/' + str(member) + '/SImon/'
else:
    dir_input = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/' + str(member) + '/SImon/'
if institute == 'CNRM-CERFACS' or model == 'MIROC-ES2L':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/r1i1p1f2/Ofx/areacello/gn/latest/'
elif model == 'CanESM5-CanOE':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/r1i1p2f1/Ofx/areacello/gn/latest/'
elif model == 'CESM2-WACCM':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/CESM2/' + str(experiment) + '/r1i1p1f1/Ofx/areacello/gn/latest/'
elif model == 'FGOALS-g3':
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/FGOALS-f3-L/' + str(experiment) + '/r1i1p1f1/Ofx/areacello/gn/latest/'
elif model == 'TaiESM1' or model == 'CAMS-CSM1-0' or institute == 'E3SM-Project' or model == 'FIO-ESM-2-0' or institute == 'INM' or institute == 'MOHC' or model == 'GISS-E2-1-H' or model == 'NESM3' or institute == 'BCC' or model == 'KIOST-ESM':
    dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/areacello/'
    if model == 'GISS-E2-1-H':
        dir_grid2 = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/areacella/'
elif model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
    dir_grid = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/areacella/'
elif model == 'CMCC-CM2-HR4':
    dir_grid = '/home/rossby/data_lib/esgf/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/r1i1p1f1/Ofx/areacello/gn/latest/'
else:
    dir_grid = '/nobackup/rossby22/sm_klazi/data_lib/CMIP6/CMIP/' + str(institute) + '/' + str(model) + '/' + str(experiment) + '/r1i1p1f1/Ofx/areacello/gn/latest/'
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
        file_grid = dir_grid + 'areacello_Ofx_' + str(model) + '_' + str(experiment) + '_r1i1p1f2_gn.nc'
    elif model == 'CanESM5-CanOE':
        file_grid = dir_grid + 'areacello_Ofx_' + str(model) + '_' + str(experiment) + '_r1i1p2f1_gn.nc'
    elif model == 'CESM2-WACCM':
        file_grid = dir_grid + 'areacello_Ofx_CESM2_' + str(experiment) + '_r1i1p1f1_gn.nc'
    elif model == 'TaiESM1' or model == 'CAMS-CSM1-0' or institute == 'E3SM-Project' or model == 'FIO-ESM-2-0' or institute == 'INM' or institute == 'MOHC' or model == 'GISS-E2-1-H' or model == 'NESM3' or model == 'KIOST-ESM':
        file_grid = dir_grid + 'areacello_' + str(model) + '.nc'
        if model == 'GISS-E2-1-H':
            file_grid2 = dir_grid2 + 'areacella_GISS-E2-1-G.nc'
    elif model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
        file_grid = dir_grid + 'areacella_GISS-E2-1-G.nc'
    elif institute == 'BCC':
        file_grid = dir_grid + 'areacello_' + str(institute) + '.nc'
    elif model == 'FGOALS-g3':
        file_grid = dir_grid + 'areacello_Ofx_FGOALS-f3-L_' + str(experiment) + '_r1i1p1f1_gn.nc'
    else:
        file_grid = dir_grid + 'areacello_Ofx_' + str(model) + '_' + str(experiment) + '_r1i1p1f1_gn.nc'
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
        if institute == 'CMCC':
            grid_area = grid_area[0:-1,1:-1]
            lat = lat[0:-1,1:-1]
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
        if (member == 'r10i1p1f1' and (start_year+year) == 1952) or (member == 'r10i1p1f1' and (start_year+year) == 1990) or (member == 'r10i1p1f1' and (start_year+year) == 2000) or (member == 'r21i1p1f1' and (start_year+year) == 1957) or (member == 'r21i1p1f1' and (start_year+year) == 1979):
            file_siconc = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/siconc/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
        else:
            file_siconc = dir_input + 'siconc/gn/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
        fh = Dataset(file_siconc, mode='r')
        siconc = fh.variables['siconc'][:]
        siconc[siconc < 0.] = np.nan
        siconc[siconc > 101.] = np.nan
        fh.close()

        # Retrieve sivol
        file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_' + str(start_year+year) + '01-' + str(start_year+year) + '12.nc'
        fh = Dataset(file_sivol, mode='r')
        sivol = fh.variables['sivol'][:]
        sivol[sivol < 0.] = np.nan
        sivol[sivol > 20.] = np.nan
        fh.close()

        # Compute total Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3)
        area_month = compute_area(nmy,lat,lat_threshold,siconc,grid_area)
        volume_month = compute_vol(nmy,lat,lat_threshold,sivol,grid_area)

        # Store variables
        for month in np.arange(nmy):
            area_total[year,month] = area_month[month]
            volume_total[year,month] = volume_month[month]

elif institute == 'HAMMOZ-Consortium' or institute == 'MPI-M' or (model == 'CESM2' and member == 'r7i1p1f1') or (model == 'CESM2' and member == 'r8i1p1f1') or (model == 'CESM2' and member == 'r9i1p1f1') or (model == 'CESM2' and member == 'r10i1p1f1') or (model == 'CESM2' and member == 'r11i1p1f1') or model == 'CESM2-FV2' or model == 'CESM2-WACCM-FV2' or model == 'NorESM2-LM' or model == 'NorESM2-MM' or model == 'SAM0-UNICON' or institute == 'E3SM-Project' or institute == 'INM' or model == 'HadGEM3-GC31-MM' or institute == 'NASA-GISS' or institute == 'AWI' or model == 'CMCC-CM2-HR4':

    # Loop over files
    if model == 'MPI-ESM1-2-HR' or institute == 'E3SM-Project':
        n_files = 13
    elif model == 'NorESM2-LM' or model == 'NorESM2-MM' or model == 'SAM0-UNICON':
        n_files = 7
    elif institute == 'NCAR' or institute == 'INM' or model == 'CMCC-CM2-HR4':
        n_files = 2
    elif model == 'HadGEM3-GC31-MM':
        n_files = 5
    elif institute == 'AWI':
        n_files = 8
    elif model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
        n_files = 3
    else:
        n_files = 4

    for filex in np.arange(n_files):
        print('file',filex+1)
        if model == 'MPI-ESM1-2-HR' or institute == 'E3SM-Project':
            nyears_file = 5
            j = nyears_file - 1
        elif model == 'NorESM2-LM' or model == 'NorESM2-MM' or model == 'SAM0-UNICON':
            nyears_file = 10
            if filex == 6:
                j = 4
            else:
                j = nyears_file - 1
        elif institute == 'NCAR' or institute == 'INM' or model == 'CMCC-CM2-HR4':
            nyears_file = 50
            if filex == 0:
                j = nyears_file - 1
            elif filex == 1:
                j = 14
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
        elif institute == 'AWI':
            if filex == 0:
                start_period = '1941'
                end_period = '1950'
                nyears_file = 1
            elif filex == 1:
                start_period = '1951'
                end_period = '1960'
                nyears_file = 10
            elif filex == 2:
                start_period = '1961'
                end_period = '1970'
                nyears_file = 10
            elif filex == 3:
                start_period = '1971'
                end_period = '1980'
                nyears_file = 10
            elif filex == 4:
                start_period = '1981'
                end_period = '1990'
                nyears_file = 10
            elif filex == 5:
                start_period = '1991'
                end_period = '2000'
                nyears_file = 10
            elif filex == 6:
                start_period = '2001'
                end_period = '2010'
                nyears_file = 10
            elif filex == 7:
                start_period = '2011'
                end_period = '2014'
                nyears_file = 4
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
            if filex == 3:
                j = 4
            else:
                j = nyears_file - 1

        # Retrieve siconc
        if not (institute == 'AWI' or model == 'GISS-E2-1-G-CC'):
            if institute == 'E3SM-Project' or model == 'GISS-E2-1-H':
                file_siconc = dir_input + 'siconc/gr/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
            elif institute == 'INM':
                file_siconc = dir_input + 'siconc/gr1/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr1_' + str(start_year+filex*nyears_file) + '01-' + str(start_year+filex*nyears_file+j) + '12.nc'
            elif model == 'HadGEM3-GC31-MM':
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
        if not (institute == 'E3SM-Project' or institute == 'INM' or institute == 'AWI' or model == 'GISS-E2-1-H'):
            if model == 'HadGEM3-GC31-MM' or model == 'GISS-E2-1-G' or model == 'GISS-E2-1-G-CC':
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
            if filex == 0:
                area_month = area_month[108::]
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
            if filex == 0:
                volume_month = volume_month[108::]
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
        if not (institute == 'INM' or model == 'GISS-E2-1-H' or institute == 'AWI'):
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
                        area_total[year+20,month] = area_month[month+year*nmy]
                        volume_total[year+20,month] = volume_month[month+year*nmy]
                    elif filex == 2:
                        area_total[year+38,month] = area_month[month+year*nmy]
                        volume_total[year+38,month] = volume_month[month+year*nmy]
                    elif filex == 3:
                        area_total[year+40,month] = area_month[month+year*nmy]
                        volume_total[year+40,month] = volume_month[month+year*nmy]
                    elif filex == 4:
                        area_total[year+60,month] = area_month[month+year*nmy]
                        volume_total[year+60,month] = volume_month[month+year*nmy]
                elif institute == 'AWI':
                    if filex == 0:
                        area_total[year,month] = area_month[month+year*nmy]
                        volume_total[year,month] = volume_month[month+year*nmy]
                    elif filex == 1:
                        area_total[year+1,month] = area_month[month+year*nmy]
                        volume_total[year+1,month] = volume_month[month+year*nmy]
                    elif filex == 2:
                        area_total[year+11,month] = area_month[month+year*nmy]
                        volume_total[year+11,month] = volume_month[month+year*nmy]
                    elif filex == 3:
                        area_total[year+21,month] = area_month[month+year*nmy]
                        volume_total[year+21,month] = volume_month[month+year*nmy]
                    elif filex == 4:
                        area_total[year+31,month] = area_month[month+year*nmy]
                        volume_total[year+31,month] = volume_month[month+year*nmy]
                    elif filex == 5:
                        area_total[year+41,month] = area_month[month+year*nmy]
                        volume_total[year+41,month] = volume_month[month+year*nmy]
                    elif filex == 6:
                        area_total[year+51,month] = area_month[month+year*nmy]
                        volume_total[year+51,month] = volume_month[month+year*nmy]
                    elif filex == 7:
                        area_total[year+61,month] = area_month[month+year*nmy]
                        volume_total[year+61,month] = volume_month[month+year*nmy]
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
                    if institute == 'INM' or model == 'GISS-E2-1-H':
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
    if model == 'MIROC6' or institute == 'NOAA-GFDL' or model == 'CIESM' or model == 'HadGEM3-GC31-LL' or model == 'UKESM1-0-LL':
        file_siconc = dir_input + 'siconc/gn/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-201412.nc'
    elif model == 'KIOST-ESM':
        file_siconc = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/historical/siconc/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr1_185001-201412.nc'
    else:
        file_siconc = dir_input + 'siconc/gn/latest/siconc_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_185001-201412.nc'
    fh = Dataset(file_siconc, mode='r')
    siconc = fh.variables['siconc'][:]
    if model == 'KIOST-ESM':
        siconc = siconc * 100.
    siconc[siconc < 0.] = np.nan
    siconc[siconc > 101.] = np.nan
    fh.close()

    # Retrieve sivol
    if not (institute == 'CCCma' or institute == 'MIROC' or model == 'NorCPM1' or model == 'FGOALS-g3' or model == 'NESM3' or model == 'KIOST-ESM'):
        if institute == 'NOAA-GFDL' or model == 'CIESM' or model == 'HadGEM3-GC31-LL' or model == 'UKESM1-0-LL':
            file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-201412.nc'
        else:
            file_sivol = dir_input + 'sivol/gn/latest/sivol_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_185001-201412.nc'
        fh = Dataset(file_sivol, mode='r')
        sivol = fh.variables['sivol'][:]
        sivol[sivol < 0.] = np.nan
        sivol[sivol > 20.] = np.nan
        fh.close()

    # Retrieve sithick and compute sivol
    if model == 'CanESM5' or institute == 'MIROC' or model == 'NESM3' or model == 'NorCPM1' or model == 'KIOST-ESM':
       if model == 'MIROC6':
           file_sithick = dir_input + 'sithick/gn/latest/sithick_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_195001-201412.nc'
       elif model == 'KIOST-ESM':
            file_sithick = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/historical/sithick/sithick_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gr1_185001-201412.nc'
       else:
           file_sithick = dir_input + 'sithick/gn/latest/sithick_SImon_' + str(model) + '_' + str(experiment) + '_' + str(member) + '_gn_185001-201412.nc'
       fh = Dataset(file_sithick, mode='r')
       sithick = fh.variables['sithick'][:]
       sithick[sithick < 0.] = np.nan
       sithick[sithick > 40.] = np.nan
       fh.close()
       sivol = (sithick * siconc) / 100.

    # Compute total Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3)
    if model == 'MIROC6' or institute == 'NOAA-GFDL' or model == 'CIESM' or model == 'HadGEM3-GC31-LL' or model == 'UKESM1-0-LL':
        area_month = compute_area(nm,lat,lat_threshold,siconc,grid_area)
    else:
        area_month = compute_area(nm,lat,lat_threshold,siconc[1200::,:,:],grid_area)
    if not (model == 'CanESM5-CanOE' or model == 'FGOALS-g3'):
        if institute == 'NOAA-GFDL' or model == 'CIESM' or model == 'HadGEM3-GC31-LL' or model == 'UKESM1-0-LL' or institute == 'MIROC':
            volume_month = compute_vol(nm,lat,lat_threshold,sivol,grid_area)
        else:
            volume_month = compute_vol(nm,lat,lat_threshold,sivol[1200::,:,:],grid_area)

    # Store variables
    for month in np.arange(nmy):
        for year in np.arange(nyears):
            area_total[year,month] = area_month[month+year*nmy]
            if model == 'CanESM5-CanOE' or model == 'FGOALS-g3':
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
