#!/usr/bin/env
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 3: Plot changes in Arctic sea-ice area (10^6 km^2) and volume (10^3 km^3) between 2015-2019 and 2096-2100 for selected CMIP6 models based on SSP5-8.5
    Fig. S5: Same as Fig. 3 for SSP1-2.6
PROGRAMMER
    D. Docquier
LAST UPDATEs
    30/04/2021
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
import random

# Option
experiment = 'ssp585' # ssp126; ssp585
save_fig = True
n_rand = 10 # number of models included in the random analysis
n_iter = 1000 # number of iterations included in the random analysis

# Time parameter
nmy = 12

# Function to compute mean SIA/SIV at beginning and end of scenario and relative change in SIA and SIV
def compute_change(var):
    start = np.nanmean(var[0:5,:],axis=0) # 5-year time period
    end = np.nanmean(var[-5::,:],axis=0)
#    start = np.nanmean(var[0:10,:],axis=0) # 10-year time period
#    end = np.nanmean(var[-10::,:],axis=0)
    change = 100. * (end - start) / start
    return start,end,change

# Function to retrieve monthly mean
def compute_monthly(nm,nmy,var):
    nyears = int(nm / nmy)
    var_mon = np.zeros((nyears,nmy))
    for k in np.arange(nyears):
        for i in np.arange(nmy):
            if np.isnan(var[i+k*nmy]) == False:
                var_mon[k,i] = var[i+k*nmy]
    return var_mon

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/' + str(experiment) + '/SIA-SIV/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6_Paper/'

# Load sea-ice area and volume AWI-CM-1-1-MR
filename = dir_input + 'SIA-SIV_SImon_AWI-CM-1-1-MR_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_awicm_init,volume_awicm_init = np.load(filename)
area_awicm_start,area_awicm_end,area_awicm_change = compute_change(area_awicm_init)
volume_awicm_start,volume_awicm_end,volume_awicm_change = compute_change(volume_awicm_init)

# Load sea-ice area and volume BCC-CSM2-MR
filename = dir_input + 'SIA-SIV_SImon_BCC-CSM2-MR_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_bcccsm2mr_init,volume_bcccsm2mr_init = np.load(filename)
area_bcccsm2mr_start,area_bcccsm2mr_end,area_bcccsm2mr_change = compute_change(area_bcccsm2mr_init)
volume_bcccsm2mr_start,volume_bcccsm2mr_end,volume_bcccsm2mr_change = compute_change(volume_bcccsm2mr_init)

# Load sea-ice area and volume CAMS-CSM1-0
filename = dir_input + 'SIA-SIV_SImon_CAMS-CSM1-0_' + str(experiment) + '_ensmean_gn_2015-2099.npy'
area_camscsm_init,volume_camscsm_init,notused,notused = np.load(filename)
area_camscsm_init2 = np.zeros((86,12)) * np.nan
volume_camscsm_init2 = np.zeros((86,12)) * np.nan
area_camscsm_init2[0:85,:] = area_camscsm_init
volume_camscsm_init2[0:85,:] = volume_camscsm_init
area_camscsm_start,area_camscsm_end,area_camscsm_change = compute_change(area_camscsm_init2)
volume_camscsm_start,volume_camscsm_end,volume_camscsm_change = compute_change(volume_camscsm_init2)

# Load sea-ice area and volume FGOALS-f3-L
filename = dir_input + 'SIA-SIV_SImon_FGOALS-f3-L_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_fgoalsf3l_init,volume_fgoalsf3l_init = np.load(filename)
area_fgoalsf3l_start,area_fgoalsf3l_end,area_fgoalsf3l_change = compute_change(area_fgoalsf3l_init)
volume_fgoalsf3l_start,volume_fgoalsf3l_end,volume_fgoalsf3l_change = compute_change(volume_fgoalsf3l_init)

# Load sea-ice area and volume FGOALS-g3
filename = dir_input + 'SIA-SIV_SImon_FGOALS-g3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_fgoalsg3_init,volume_fgoalsg3_init,notused,notused = np.load(filename)
area_fgoalsg3_start,area_fgoalsg3_end,area_fgoalsg3_change = compute_change(area_fgoalsg3_init)
volume_fgoalsg3_start,volume_fgoalsg3_end,volume_fgoalsg3_change = compute_change(volume_fgoalsg3_init)

# Load sea-ice area and volume CanESM5
filename = dir_input + 'SIA-SIV_SImon_CanESM5_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_canesm5_init,volume_canesm5_init,notused,notused = np.load(filename)
area_canesm5_start,area_canesm5_end,area_canesm5_change = compute_change(area_canesm5_init)
volume_canesm5_start,volume_canesm5_end,volume_canesm5_change = compute_change(volume_canesm5_init)

# Load sea-ice area and volume CanESM5-CanOE
filename = dir_input + 'SIA-SIV_SImon_CanESM5-CanOE_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_canesm5canoe_init,volume_canesm5canoe_init,notused,notused = np.load(filename)
area_canesm5canoe_start,area_canesm5canoe_end,area_canesm5canoe_change = compute_change(area_canesm5canoe_init)
volume_canesm5canoe_start,volume_canesm5canoe_end,volume_canesm5canoe_change = compute_change(volume_canesm5canoe_init)

# Load sea-ice area and volume CNRM-CM6-1
filename = dir_input + 'SIA-SIV_SImon_CNRM-CM6-1_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cnrmcm6_init,volume_cnrmcm6_init,notused,notused = np.load(filename)
area_cnrmcm6_start,area_cnrmcm6_end,area_cnrmcm6_change = compute_change(area_cnrmcm6_init)
volume_cnrmcm6_start,volume_cnrmcm6_end,volume_cnrmcm6_change = compute_change(volume_cnrmcm6_init)

# Load sea-ice area and volume CNRM-CM6-1-HR
filename = dir_input + 'SIA-SIV_SImon_CNRM-CM6-1-HR_' + str(experiment) + '_r1i1p1f2_gn_2015-2100.npy'
area_cnrmcm6hr_init,volume_cnrmcm6hr_init = np.load(filename)
area_cnrmcm6hr_start,area_cnrmcm6hr_end,area_cnrmcm6hr_change = compute_change(area_cnrmcm6hr_init)
volume_cnrmcm6hr_start,volume_cnrmcm6hr_end,volume_cnrmcm6hr_change = compute_change(volume_cnrmcm6hr_init)

# Load sea-ice area and volume ACCESS-ESM1-5
filename = dir_input + 'SIA-SIV_SImon_ACCESS-ESM1-5_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_accessesm_init,volume_accessesm_init,notused,notused = np.load(filename)
area_accessesm_start,area_accessesm_end,area_accessesm_change = compute_change(area_accessesm_init)
volume_accessesm_start,volume_accessesm_end,volume_accessesm_change = compute_change(volume_accessesm_init)

# Load sea-ice area and volume ACCESS-CM2
filename = dir_input + 'SIA-SIV_SImon_ACCESS-CM2_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_accesscm_init,volume_accesscm_init,notused,notused = np.load(filename)
area_accesscm_start,area_accesscm_end,area_accesscm_change = compute_change(area_accesscm_init)
volume_accesscm_start,volume_accesscm_end,volume_accesscm_change = compute_change(volume_accesscm_init)

# Load sea-ice area and volume MPI-ESM1-2-LR
filename = dir_input + 'SIA-SIV_SImon_MPI-ESM1-2-LR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_mpiesmlr_init,volume_mpiesmlr_init,notused,notused = np.load(filename)
area_mpiesmlr_start,area_mpiesmlr_end,area_mpiesmlr_change = compute_change(area_mpiesmlr_init)
volume_mpiesmlr_start,volume_mpiesmlr_end,volume_mpiesmlr_change = compute_change(volume_mpiesmlr_init)

# Load sea-ice area and volume MPI-ESM1-2-HR
filename = dir_input + 'SIA-SIV_SImon_MPI-ESM1-2-HR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_mpiesmhr_init,volume_mpiesmhr_init,notused,notused = np.load(filename)
area_mpiesmhr_start,area_mpiesmhr_end,area_mpiesmhr_change = compute_change(area_mpiesmhr_init)
volume_mpiesmhr_start,volume_mpiesmhr_end,volume_mpiesmhr_change = compute_change(volume_mpiesmhr_init)

# Load sea-ice area and volume EC-Earth3
filename = dir_input + 'SIA-SIV_SImon_EC-Earth3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ecearth3_init,volume_ecearth3_init,notused,notused = np.load(filename)
area_ecearth3_start,area_ecearth3_end,area_ecearth3_change = compute_change(area_ecearth3_init)
volume_ecearth3_start,volume_ecearth3_end,volume_ecearth3_change = compute_change(volume_ecearth3_init)

# Load sea-ice area and volume EC-Earth3-Veg
filename = dir_input + 'SIA-SIV_SImon_EC-Earth3-Veg_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ecearth3veg_init,volume_ecearth3veg_init,notused,notused = np.load(filename)
area_ecearth3veg_start,area_ecearth3veg_end,area_ecearth3veg_change = compute_change(area_ecearth3veg_init)
volume_ecearth3veg_start,volume_ecearth3veg_end,volume_ecearth3veg_change = compute_change(volume_ecearth3veg_init)

# Load sea-ice area and volume FIO-ESM-2-0
filename = dir_input + 'SIA-SIV_SImon_FIO-ESM-2-0_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_fioesm_init,volume_fioesm_init,notused,notused = np.load(filename)
area_fioesm_start,area_fioesm_end,area_fioesm_change = compute_change(area_fioesm_init)
volume_fioesm_start,volume_fioesm_end,volume_fioesm_change = compute_change(volume_fioesm_init)

# Load sea-ice area and volume INM-CM4-8
filename = dir_input + 'SIA-SIV_SImon_INM-CM4-8_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_inmcm48_init,volume_inmcm48_init = np.load(filename)
area_inmcm48_start,area_inmcm48_end,area_inmcm48_change = compute_change(area_inmcm48_init)
volume_inmcm48_start,volume_inmcm48_end,volume_inmcm48_change = compute_change(volume_inmcm48_init)

# Load sea-ice area and volume INM-CM5-0
filename = dir_input + 'SIA-SIV_SImon_INM-CM5-0_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_inmcm50_init,volume_inmcm50_init = np.load(filename)
area_inmcm50_start,area_inmcm50_end,area_inmcm50_change = compute_change(area_inmcm50_init)
volume_inmcm50_start,volume_inmcm50_end,volume_inmcm50_change = compute_change(volume_inmcm50_init)

# Load sea-ice area and volume IPSL-CM6A-LR
filename = dir_input + 'SIA-SIV_SImon_IPSL-CM6A-LR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ipslcm6alr_init,volume_ipslcm6alr_init,notused,notused = np.load(filename)
area_ipslcm6alr_start,area_ipslcm6alr_end,area_ipslcm6alr_change = compute_change(area_ipslcm6alr_init)
volume_ipslcm6alr_start,volume_ipslcm6alr_end,volume_ipslcm6alr_change = compute_change(volume_ipslcm6alr_init)

# Load sea-ice area and volume KIOST-ESM
filename = dir_input + 'SIA-SIV_SImon_KIOST-ESM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_kiostesm_init,volume_kiostesm_init = np.load(filename)
area_kiostesm_start,area_kiostesm_end,area_kiostesm_change = compute_change(area_kiostesm_init)
volume_kiostesm_start,volume_kiostesm_end,volume_kiostesm_change = compute_change(volume_kiostesm_init)

# Load sea-ice area and volume MIROC6
filename = dir_input + 'SIA-SIV_SImon_MIROC6_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_miroc6_init,volume_miroc6_init,notused,notused = np.load(filename)
area_miroc6_start,area_miroc6_end,area_miroc6_change = compute_change(area_miroc6_init)
volume_miroc6_start,volume_miroc6_end,volume_miroc6_change = compute_change(volume_miroc6_init)

# Load sea-ice area and volume MIROC-ES2L
filename = dir_input + 'SIA-SIV_SImon_MIROC-ES2L_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_miroces2l_init,volume_miroces2l_init,notused,notused = np.load(filename)
area_miroces2l_start,area_miroces2l_end,area_miroces2l_change = compute_change(area_miroces2l_init)
volume_miroces2l_start,volume_miroces2l_end,volume_miroces2l_change = compute_change(volume_miroces2l_init)

# Load sea-ice area and volume HadGEM3-GC31-LL
filename = dir_input + 'SIA-SIV_SImon_HadGEM3-GC31-LL_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_hadgem3ll_init,volume_hadgem3ll_init,notused,notused = np.load(filename)
area_hadgem3ll_start,area_hadgem3ll_end,area_hadgem3ll_change = compute_change(area_hadgem3ll_init)
volume_hadgem3ll_start,volume_hadgem3ll_end,volume_hadgem3ll_change = compute_change(volume_hadgem3ll_init)

# Load sea-ice area and volume HadGEM3-GC31-MM
filename = dir_input + 'SIA-SIV_SImon_HadGEM3-GC31-MM_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_hadgem3mm_init,volume_hadgem3mm_init,notused,notused = np.load(filename)
area_hadgem3mm_start,area_hadgem3mm_end,area_hadgem3mm_change = compute_change(area_hadgem3mm_init)
volume_hadgem3mm_start,volume_hadgem3mm_end,volume_hadgem3mm_change = compute_change(volume_hadgem3mm_init)

# Load sea-ice area and volume UKESM1-0-LL
filename = dir_input + 'SIA-SIV_SImon_UKESM1-0-LL_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ukesmll_init,volume_ukesmll_init,notused,notused = np.load(filename)
area_ukesmll_start,area_ukesmll_end,area_ukesmll_change = compute_change(area_ukesmll_init)
volume_ukesmll_start,volume_ukesmll_end,volume_ukesmll_change = compute_change(volume_ukesmll_init)

# Load sea-ice area and volume MRI-ESM2-0
filename = dir_input + 'SIA-SIV_SImon_MRI-ESM2-0_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_mriesm_init,volume_mriesm_init = np.load(filename)
area_mriesm_start,area_mriesm_end,area_mriesm_change = compute_change(area_mriesm_init)
volume_mriesm_start,volume_mriesm_end,volume_mriesm_change = compute_change(volume_mriesm_init)

# Load sea-ice area and volume CESM2
filename = dir_input + 'SIA-SIV_SImon_CESM2_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cesm2_init,volume_cesm2_init,notused,notused = np.load(filename)
area_cesm2_start,area_cesm2_end,area_cesm2_change = compute_change(area_cesm2_init)
volume_cesm2_start,volume_cesm2_end,volume_cesm2_change = compute_change(volume_cesm2_init)

# Load sea-ice area and volume CESM2-WACCM
filename = dir_input + 'SIA-SIV_SImon_CESM2-WACCM_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cesm2waccm_init,volume_cesm2waccm_init,notused,notused = np.load(filename)
area_cesm2waccm_start,area_cesm2waccm_end,area_cesm2waccm_change = compute_change(area_cesm2waccm_init)
volume_cesm2waccm_start,volume_cesm2waccm_end,volume_cesm2waccm_change = compute_change(volume_cesm2waccm_init)

# Load sea-ice area and volume NorESM2-LM
filename = dir_input + 'SIA-SIV_SImon_NorESM2-LM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_noresm2lm_init,volume_noresm2lm_init = np.load(filename)
area_noresm2lm_start,area_noresm2lm_end,area_noresm2lm_change = compute_change(area_noresm2lm_init)
volume_noresm2lm_start,volume_noresm2lm_end,volume_noresm2lm_change = compute_change(volume_noresm2lm_init)

# Load sea-ice area and volume NorESM2-MM
filename = dir_input + 'SIA-SIV_SImon_NorESM2-MM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_noresm2mm_init,volume_noresm2mm_init = np.load(filename)
area_noresm2mm_start,area_noresm2mm_end,area_noresm2mm_change = compute_change(area_noresm2mm_init)
volume_noresm2mm_start,volume_noresm2mm_end,volume_noresm2mm_change = compute_change(volume_noresm2mm_init)

# Load sea-ice area and volume GFDL-CM4
if experiment == 'ssp585':
    filename = dir_input + 'SIA-SIV_SImon_GFDL-CM4_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
    area_gfdlcm4_init,volume_gfdlcm4_init = np.load(filename)
    area_gfdlcm4_start,area_gfdlcm4_end,area_gfdlcm4_change = compute_change(area_gfdlcm4_init)
    volume_gfdlcm4_start,volume_gfdlcm4_end,volume_gfdlcm4_change = compute_change(volume_gfdlcm4_init)

# Load sea-ice area and volume GFDL-ESM4
filename = dir_input + 'SIA-SIV_SImon_GFDL-ESM4_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_gfdlesm4_init,volume_gfdlesm4_init = np.load(filename)
area_gfdlesm4_start,area_gfdlesm4_end,area_gfdlesm4_change = compute_change(area_gfdlesm4_init)
volume_gfdlesm4_start,volume_gfdlesm4_end,volume_gfdlesm4_change = compute_change(volume_gfdlesm4_init)

# Load sea-ice area and volume NESM3
filename = dir_input + 'SIA-SIV_SImon_NESM3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_nesm3_init,volume_nesm3_init,notused,notused = np.load(filename)
area_nesm3_start,area_nesm3_end,area_nesm3_change = compute_change(area_nesm3_init)
volume_nesm3_start,volume_nesm3_end,volume_nesm3_change = compute_change(volume_nesm3_init)

# Load Arctic sea-ice area OSI SAF (OSI-450) 1979-2015
filename = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/post-proc/t613/SIarea_OSI-450_1979-2015.npy'
area_obs1,notused,notused,notused,notused,notused,notused,notused = np.load(filename)

# Load Arctic sea-ice area OSI SAF (OSI-430-b) from 2016
filename = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/obs/OSI-430-b/SIarea_OSI-430-b_2016-2020.npy'
area_obs2,notused,notused,notused,notused,notused,notused,notused,notused = np.load(filename)

# Merge observed area 1979-2020, separate by month
area_obs_init = np.concatenate((area_obs1,area_obs2))
nm = np.size(area_obs_init)
area_obs = compute_monthly(nm,nmy,area_obs_init)
area_obs_start = np.nanmean(area_obs[-6:-1,:],axis=0)

# Load sea-ice volume PIOMAS 1979-2019
year,vol1,vol2,vol3,vol4,vol5,vol6,vol7,vol8,vol9,vol10,vol11,vol12 = np.loadtxt('/nobackup/rossby24/proj/rossby/joint_exp/oseaice/obs/PIOMAS.2sst.monthly.Current.v2.1.txt',unpack=True)
nyears_piomas = np.size(year) - 1 # Remove 2020, which is incomplete
nm_piomas = nyears_piomas * nmy
volume_piomas = np.zeros((nyears_piomas,nmy))
k = 0
for year in np.arange(nyears_piomas):
    volume_piomas[year,k] = vol1[year]
    volume_piomas[year,k+1] = vol2[year]
    volume_piomas[year,k+2] = vol3[year]
    volume_piomas[year,k+3] = vol4[year]
    volume_piomas[year,k+4] = vol5[year]
    volume_piomas[year,k+5] = vol6[year]
    volume_piomas[year,k+6] = vol7[year]
    volume_piomas[year,k+7] = vol8[year]
    volume_piomas[year,k+8] = vol9[year]
    volume_piomas[year,k+9] = vol10[year]
    volume_piomas[year,k+10] = vol11[year]
    volume_piomas[year,k+11] = vol12[year]
volume_piomas_start = np.nanmean(volume_piomas[-5::,:],axis=0)

# Compute multi-model mean - all models - Start
area_mmm_start = np.zeros(nmy)
volume_mmm_start = np.zeros(nmy)
sd_area_mmm_start = np.zeros(nmy)
sd_volume_mmm_start = np.zeros(nmy)
area_select_random_start = np.zeros(nmy)
volume_select_random_start = np.zeros(nmy)
sd_area_select_random_start = np.zeros(nmy)
sd_volume_select_random_start = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_awicm_start[m],area_bcccsm2mr_start[m],area_camscsm_start[m],area_fgoalsf3l_start[m],area_fgoalsg3_start[m],area_canesm5_start[m],area_canesm5canoe_start[m],area_cnrmcm6_start[m],area_cnrmcm6hr_start[m],area_accessesm_start[m],area_accesscm_start[m],area_mpiesmlr_start[m],area_mpiesmhr_start[m],area_ecearth3_start[m],area_ecearth3veg_start[m],area_fioesm_start[m],area_inmcm48_start[m],area_inmcm50_start[m],area_ipslcm6alr_start[m],area_kiostesm_start[m],area_miroc6_start[m],area_miroces2l_start[m],area_hadgem3ll_start[m],area_hadgem3mm_start[m],area_ukesmll_start[m],area_mriesm_start[m],area_cesm2_start[m],area_cesm2waccm_start[m],area_noresm2lm_start[m],area_noresm2mm_start[m],area_gfdlcm4_start[m],area_gfdlesm4_start[m],area_nesm3_start[m]]
        array_volume = [volume_awicm_start[m],volume_bcccsm2mr_start[m],volume_camscsm_start[m],volume_fgoalsf3l_start[m],volume_canesm5_start[m],volume_cnrmcm6_start[m],volume_cnrmcm6hr_start[m],volume_accessesm_start[m],volume_accesscm_start[m],volume_mpiesmlr_start[m],volume_mpiesmhr_start[m],volume_ecearth3_start[m],volume_ecearth3veg_start[m],volume_fioesm_start[m],volume_ipslcm6alr_start[m],volume_kiostesm_start[m],volume_miroc6_start[m],volume_miroces2l_start[m],volume_hadgem3ll_start[m],volume_hadgem3mm_start[m],volume_ukesmll_start[m],volume_mriesm_start[m],volume_cesm2_start[m],volume_cesm2waccm_start[m],volume_noresm2lm_start[m],volume_noresm2mm_start[m],volume_gfdlcm4_start[m],volume_nesm3_start[m]]
    else:
        array_area = [area_awicm_start[m],area_bcccsm2mr_start[m],area_camscsm_start[m],area_fgoalsf3l_start[m],area_fgoalsg3_start[m],area_canesm5_start[m],area_canesm5canoe_start[m],area_cnrmcm6_start[m],area_cnrmcm6hr_start[m],area_accessesm_start[m],area_accesscm_start[m],area_mpiesmlr_start[m],area_mpiesmhr_start[m],area_ecearth3_start[m],area_ecearth3veg_start[m],area_fioesm_start[m],area_inmcm48_start[m],area_inmcm50_start[m],area_ipslcm6alr_start[m],area_kiostesm_start[m],area_miroc6_start[m],area_miroces2l_start[m],area_hadgem3ll_start[m],area_hadgem3mm_start[m],area_ukesmll_start[m],area_mriesm_start[m],area_cesm2_start[m],area_cesm2waccm_start[m],area_noresm2lm_start[m],area_noresm2mm_start[m],area_gfdlesm4_start[m],area_nesm3_start[m]]
        array_volume = [volume_awicm_start[m],volume_bcccsm2mr_start[m],volume_camscsm_start[m],volume_fgoalsf3l_start[m],volume_canesm5_start[m],volume_cnrmcm6_start[m],volume_cnrmcm6hr_start[m],volume_accessesm_start[m],volume_accesscm_start[m],volume_mpiesmlr_start[m],volume_mpiesmhr_start[m],volume_ecearth3_start[m],volume_ecearth3veg_start[m],volume_fioesm_start[m],volume_ipslcm6alr_start[m],volume_kiostesm_start[m],volume_miroc6_start[m],volume_miroces2l_start[m],volume_hadgem3ll_start[m],volume_hadgem3mm_start[m],volume_ukesmll_start[m],volume_mriesm_start[m],volume_cesm2_start[m],volume_cesm2waccm_start[m],volume_noresm2lm_start[m],volume_noresm2mm_start[m],volume_gfdlesm4_start[m],volume_nesm3_start[m]]
    area_mmm_start[m] = np.nanmean(array_area)
    volume_mmm_start[m] = np.nanmean(array_volume)
    sd_area_mmm_start[m] = np.nanstd(array_area)
    sd_volume_mmm_start[m] = np.nanstd(array_volume)

    # Make random analysis
    area_random_mean = np.zeros(n_iter)
    volume_random_mean = np.zeros(n_iter)
    for i in np.arange(n_iter):
        array_area_random = random.sample(array_area,n_rand)
        array_volume_random = random.sample(array_volume,n_rand)
        area_random_mean[i] = np.nanmean(array_area_random)
        volume_random_mean[i] = np.nanmean(array_volume_random)
    area_select_random_start[m] = np.nanmean(area_random_mean)
    volume_select_random_start[m] = np.nanmean(volume_random_mean)
    sd_area_select_random_start[m] = np.nanstd(area_random_mean)
    sd_volume_select_random_start[m] = np.nanstd(volume_random_mean)

# Compute multi-model mean - all models - End
area_mmm_end = np.zeros(nmy)
volume_mmm_end = np.zeros(nmy)
sd_area_mmm_end = np.zeros(nmy)
sd_volume_mmm_end = np.zeros(nmy)
area_mmm_change = np.zeros(nmy)
volume_mmm_change = np.zeros(nmy)
area_select_random_end = np.zeros(nmy)
volume_select_random_end = np.zeros(nmy)
sd_area_select_random_end = np.zeros(nmy)
sd_volume_select_random_end = np.zeros(nmy)
area_select_random_change = np.zeros(nmy)
volume_select_random_change = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_awicm_end[m],area_bcccsm2mr_end[m],area_camscsm_end[m],area_fgoalsf3l_end[m],area_fgoalsg3_end[m],area_canesm5_end[m],area_canesm5canoe_end[m],area_cnrmcm6_end[m],area_cnrmcm6hr_end[m],area_accessesm_end[m],area_accesscm_end[m],area_mpiesmlr_end[m],area_mpiesmhr_end[m],area_ecearth3_end[m],area_ecearth3veg_end[m],area_fioesm_end[m],area_inmcm48_end[m],area_inmcm50_end[m],area_ipslcm6alr_end[m],area_kiostesm_end[m],area_miroc6_end[m],area_miroces2l_end[m],area_hadgem3ll_end[m],area_hadgem3mm_end[m],area_ukesmll_end[m],area_mriesm_end[m],area_cesm2_end[m],area_cesm2waccm_end[m],area_noresm2lm_end[m],area_noresm2mm_end[m],area_gfdlcm4_end[m],area_gfdlesm4_end[m],area_nesm3_end[m]]
        array_volume = [volume_awicm_end[m],volume_bcccsm2mr_end[m],volume_camscsm_end[m],volume_fgoalsf3l_end[m],volume_canesm5_end[m],volume_cnrmcm6_end[m],volume_cnrmcm6hr_end[m],volume_accessesm_end[m],volume_accesscm_end[m],volume_mpiesmlr_end[m],volume_mpiesmhr_end[m],volume_ecearth3_end[m],volume_ecearth3veg_end[m],volume_fioesm_end[m],volume_ipslcm6alr_end[m],volume_kiostesm_end[m],volume_miroc6_end[m],volume_miroces2l_end[m],volume_hadgem3ll_end[m],volume_hadgem3mm_end[m],volume_ukesmll_end[m],volume_mriesm_end[m],volume_cesm2_end[m],volume_cesm2waccm_end[m],volume_noresm2lm_end[m],volume_noresm2mm_end[m],volume_gfdlcm4_end[m],volume_nesm3_end[m]]
    else:
        array_area = [area_awicm_end[m],area_bcccsm2mr_end[m],area_camscsm_end[m],area_fgoalsf3l_end[m],area_fgoalsg3_end[m],area_canesm5_end[m],area_canesm5canoe_end[m],area_cnrmcm6_end[m],area_cnrmcm6hr_end[m],area_accessesm_end[m],area_accesscm_end[m],area_mpiesmlr_end[m],area_mpiesmhr_end[m],area_ecearth3_end[m],area_ecearth3veg_end[m],area_fioesm_end[m],area_inmcm48_end[m],area_inmcm50_end[m],area_ipslcm6alr_end[m],area_kiostesm_end[m],area_miroc6_end[m],area_miroces2l_end[m],area_hadgem3ll_end[m],area_hadgem3mm_end[m],area_ukesmll_end[m],area_mriesm_end[m],area_cesm2_end[m],area_cesm2waccm_end[m],area_noresm2lm_end[m],area_noresm2mm_end[m],area_gfdlesm4_end[m],area_nesm3_end[m]]
        array_volume = [volume_awicm_end[m],volume_bcccsm2mr_end[m],volume_camscsm_end[m],volume_fgoalsf3l_end[m],volume_canesm5_end[m],volume_cnrmcm6_end[m],volume_cnrmcm6hr_end[m],volume_accessesm_end[m],volume_accesscm_end[m],volume_mpiesmlr_end[m],volume_mpiesmhr_end[m],volume_ecearth3_end[m],volume_ecearth3veg_end[m],volume_fioesm_end[m],volume_ipslcm6alr_end[m],volume_kiostesm_end[m],volume_miroc6_end[m],volume_miroces2l_end[m],volume_hadgem3ll_end[m],volume_hadgem3mm_end[m],volume_ukesmll_end[m],volume_mriesm_end[m],volume_cesm2_end[m],volume_cesm2waccm_end[m],volume_noresm2lm_end[m],volume_noresm2mm_end[m],volume_gfdlesm4_end[m],volume_nesm3_end[m]]
    area_mmm_end[m] = np.nanmean(array_area)
    volume_mmm_end[m] = np.nanmean(array_volume)
    sd_area_mmm_end[m] = np.nanstd(array_area)
    sd_volume_mmm_end[m] = np.nanstd(array_volume)
    area_mmm_change[m] = 100. * (area_mmm_end[m] - area_mmm_start[m]) / area_mmm_start[m]
    volume_mmm_change[m] = 100. * (volume_mmm_end[m] - volume_mmm_start[m]) / volume_mmm_start[m]
    
    # Make random analysis
    area_random_mean = np.zeros(n_iter)
    volume_random_mean = np.zeros(n_iter)
    for i in np.arange(n_iter):
        array_area_random = random.sample(array_area,n_rand)
        array_volume_random = random.sample(array_volume,n_rand)
        area_random_mean[i] = np.nanmean(array_area_random)
        volume_random_mean[i] = np.nanmean(array_volume_random)
    area_select_random_end[m] = np.nanmean(area_random_mean)
    volume_select_random_end[m] = np.nanmean(volume_random_mean)
    sd_area_select_random_end[m] = np.nanstd(area_random_mean)
    sd_volume_select_random_end[m] = np.nanstd(volume_random_mean)
    area_select_random_change[m] = 100. * (area_select_random_end[m] - area_select_random_start[m]) / area_select_random_start[m]
    volume_select_random_change[m] = 100. * (volume_select_random_end[m] - volume_select_random_start[m]) / volume_select_random_start[m]

# Compute multi-model mean - good mean SIA (15 best models) - Start
area_select_sia_15_start = np.zeros(nmy)
volume_select_sia_15_start = np.zeros(nmy)
sd_area_select_sia_15_start = np.zeros(nmy)
sd_volume_select_sia_15_start = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_cesm2waccm_start[m],area_noresm2lm_start[m],area_accessesm_start[m],area_gfdlesm4_start[m],area_mpiesmlr_start[m],area_cnrmcm6hr_start[m],area_hadgem3mm_start[m],area_ipslcm6alr_start[m],area_gfdlcm4_start[m],area_inmcm50_start[m],area_cnrmcm6_start[m],area_ecearth3_start[m],area_ecearth3veg_start[m],area_hadgem3ll_start[m],area_miroces2l_start[m]]
        array_volume = [volume_cesm2waccm_start[m],volume_noresm2lm_start[m],volume_accessesm_start[m],volume_gfdlesm4_start[m],volume_mpiesmlr_start[m],volume_cnrmcm6hr_start[m],volume_hadgem3mm_start[m],volume_ipslcm6alr_start[m],volume_gfdlcm4_start[m],volume_cnrmcm6_start[m],volume_ecearth3_start[m],volume_ecearth3veg_start[m],volume_hadgem3ll_start[m],volume_miroces2l_start[m],volume_mriesm_start[m]]
    else:
        array_area = [area_cesm2waccm_start[m],area_noresm2lm_start[m],area_accessesm_start[m],area_gfdlesm4_start[m],area_mpiesmlr_start[m],area_cnrmcm6hr_start[m],area_hadgem3mm_start[m],area_ipslcm6alr_start[m],area_inmcm50_start[m],area_cnrmcm6_start[m],area_ecearth3_start[m],area_ecearth3veg_start[m],area_hadgem3ll_start[m],area_miroces2l_start[m],area_mriesm_start[m]]
        array_volume = [volume_cesm2waccm_start[m],volume_noresm2lm_start[m],volume_accessesm_start[m],volume_gfdlesm4_start[m],volume_mpiesmlr_start[m],volume_cnrmcm6hr_start[m],volume_hadgem3mm_start[m],volume_ipslcm6alr_start[m],volume_cnrmcm6_start[m],volume_ecearth3_start[m],volume_ecearth3veg_start[m],volume_hadgem3ll_start[m],volume_miroces2l_start[m],volume_mriesm_start[m],volume_accesscm_start[m]]
    area_select_sia_15_start[m] = np.nanmean(array_area)
    volume_select_sia_15_start[m] = np.nanmean(array_volume)
    sd_area_select_sia_15_start[m] = np.nanstd(array_area)
    sd_volume_select_sia_15_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good mean SIA (15 best models) - End
area_select_sia_15_end = np.zeros(nmy)
volume_select_sia_15_end = np.zeros(nmy)
sd_area_select_sia_15_end = np.zeros(nmy)
sd_volume_select_sia_15_end = np.zeros(nmy)
area_select_sia_15_change = np.zeros(nmy)
volume_select_sia_15_change = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_cesm2waccm_end[m],area_noresm2lm_end[m],area_accessesm_end[m],area_gfdlesm4_end[m],area_mpiesmlr_end[m],area_cnrmcm6hr_end[m],area_hadgem3mm_end[m],area_ipslcm6alr_end[m],area_gfdlcm4_end[m],area_inmcm50_end[m],area_cnrmcm6_end[m],area_ecearth3_end[m],area_ecearth3veg_end[m],area_hadgem3ll_end[m],area_miroces2l_end[m]]
        array_volume = [volume_cesm2waccm_end[m],volume_noresm2lm_end[m],volume_accessesm_end[m],volume_gfdlesm4_end[m],volume_mpiesmlr_end[m],volume_cnrmcm6hr_end[m],volume_hadgem3mm_end[m],volume_ipslcm6alr_end[m],volume_gfdlcm4_end[m],volume_cnrmcm6_end[m],volume_ecearth3_end[m],volume_ecearth3veg_end[m],volume_hadgem3ll_end[m],volume_miroces2l_end[m],volume_mriesm_end[m]]
    else:
        array_area = [area_cesm2waccm_end[m],area_noresm2lm_end[m],area_accessesm_end[m],area_gfdlesm4_end[m],area_mpiesmlr_end[m],area_cnrmcm6hr_end[m],area_hadgem3mm_end[m],area_ipslcm6alr_end[m],area_inmcm50_end[m],area_cnrmcm6_end[m],area_ecearth3_end[m],area_ecearth3veg_end[m],area_hadgem3ll_end[m],area_miroces2l_end[m],area_mriesm_end[m]]
        array_volume = [volume_cesm2waccm_end[m],volume_noresm2lm_end[m],volume_accessesm_end[m],volume_gfdlesm4_end[m],volume_mpiesmlr_end[m],volume_cnrmcm6hr_end[m],volume_hadgem3mm_end[m],volume_ipslcm6alr_end[m],volume_cnrmcm6_end[m],volume_ecearth3_end[m],volume_ecearth3veg_end[m],volume_hadgem3ll_end[m],volume_miroces2l_end[m],volume_mriesm_end[m],volume_accesscm_end[m]]
    area_select_sia_15_end[m] = np.nanmean(array_area)
    volume_select_sia_15_end[m] = np.nanmean(array_volume)
    sd_area_select_sia_15_end[m] = np.nanstd(array_area)
    sd_volume_select_sia_15_end[m] = np.nanstd(array_volume)
    area_select_sia_15_change[m] = 100. * (area_select_sia_15_end[m] - area_select_sia_15_start[m]) / area_select_sia_15_start[m]
    volume_select_sia_15_change[m] = 100. * (volume_select_sia_15_end[m] - volume_select_sia_15_start[m]) / volume_select_sia_15_start[m]

# Compute multi-model mean - good mean SIV (15 best models) - Start
area_select_siv_15_start = np.zeros(nmy)
volume_select_siv_15_start = np.zeros(nmy)
sd_area_select_siv_15_start = np.zeros(nmy)
sd_volume_select_siv_15_start = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_cesm2waccm_start[m],area_ipslcm6alr_start[m],area_fgoalsf3l_start[m],area_camscsm_start[m],area_miroc6_start[m],area_mpiesmhr_start[m],area_fioesm_start[m],area_accessesm_start[m],area_mpiesmlr_start[m],area_awicm_start[m],area_cesm2_start[m],area_gfdlcm4_start[m],area_hadgem3mm_start[m],area_mriesm_start[m],area_accesscm_start[m]]
        array_volume = [volume_cesm2waccm_start[m],volume_ipslcm6alr_start[m],volume_fgoalsf3l_start[m],volume_camscsm_start[m],volume_miroc6_start[m],volume_mpiesmhr_start[m],volume_fioesm_start[m],volume_accessesm_start[m],volume_mpiesmlr_start[m],volume_awicm_start[m],volume_cesm2_start[m],volume_gfdlcm4_start[m],volume_hadgem3mm_start[m],volume_mriesm_start[m],volume_accesscm_start[m]]
    else:
        array_area = [area_cesm2waccm_start[m],area_ipslcm6alr_start[m],area_fgoalsf3l_start[m],area_camscsm_start[m],area_miroc6_start[m],area_mpiesmhr_start[m],area_fioesm_start[m],area_accessesm_start[m],area_mpiesmlr_start[m],area_awicm_start[m],area_cesm2_start[m],area_hadgem3mm_start[m],area_mriesm_start[m],area_accesscm_start[m],area_noresm2lm_start[m]]
        array_volume = [volume_cesm2waccm_start[m],volume_ipslcm6alr_start[m],volume_fgoalsf3l_start[m],volume_camscsm_start[m],volume_miroc6_start[m],volume_mpiesmhr_start[m],volume_fioesm_start[m],volume_accessesm_start[m],volume_mpiesmlr_start[m],volume_awicm_start[m],volume_cesm2_start[m],volume_hadgem3mm_start[m],volume_mriesm_start[m],volume_accesscm_start[m],volume_noresm2lm_start[m]]
    area_select_siv_15_start[m] = np.nanmean(array_area)
    volume_select_siv_15_start[m] = np.nanmean(array_volume)
    sd_area_select_siv_15_start[m] = np.nanstd(array_area)
    sd_volume_select_siv_15_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good mean SIV (15 best models) - End
area_select_siv_15_end = np.zeros(nmy)
volume_select_siv_15_end = np.zeros(nmy)
sd_area_select_siv_15_end = np.zeros(nmy)
sd_volume_select_siv_15_end = np.zeros(nmy)
area_select_siv_15_change = np.zeros(nmy)
volume_select_siv_15_change = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_cesm2waccm_end[m],area_ipslcm6alr_end[m],area_fgoalsf3l_end[m],area_camscsm_end[m],area_miroc6_end[m],area_mpiesmhr_end[m],area_fioesm_end[m],area_accessesm_end[m],area_mpiesmlr_end[m],area_awicm_end[m],area_cesm2_end[m],area_gfdlcm4_end[m],area_hadgem3mm_end[m],area_mriesm_end[m],area_accesscm_end[m]]
        array_volume = [volume_cesm2waccm_end[m],volume_ipslcm6alr_end[m],volume_fgoalsf3l_end[m],volume_camscsm_end[m],volume_miroc6_end[m],volume_mpiesmhr_end[m],volume_fioesm_end[m],volume_accessesm_end[m],volume_mpiesmlr_end[m],volume_awicm_end[m],volume_cesm2_end[m],volume_gfdlcm4_end[m],volume_hadgem3mm_end[m],volume_mriesm_end[m],volume_accesscm_end[m]]
    else:
        array_area = [area_cesm2waccm_end[m],area_ipslcm6alr_end[m],area_fgoalsf3l_end[m],area_camscsm_end[m],area_miroc6_end[m],area_mpiesmhr_end[m],area_fioesm_end[m],area_accessesm_end[m],area_mpiesmlr_end[m],area_awicm_end[m],area_cesm2_end[m],area_hadgem3mm_end[m],area_mriesm_end[m],area_accesscm_end[m],area_noresm2lm_end[m]]
        array_volume = [volume_cesm2waccm_end[m],volume_ipslcm6alr_end[m],volume_fgoalsf3l_end[m],volume_camscsm_end[m],volume_miroc6_end[m],volume_mpiesmhr_end[m],volume_fioesm_end[m],volume_accessesm_end[m],volume_mpiesmlr_end[m],volume_awicm_end[m],volume_cesm2_end[m],volume_hadgem3mm_end[m],volume_mriesm_end[m],volume_accesscm_end[m],volume_noresm2lm_end[m]]
    area_select_siv_15_end[m] = np.nanmean(array_area)
    volume_select_siv_15_end[m] = np.nanmean(array_volume)
    sd_area_select_siv_15_end[m] = np.nanstd(array_area)
    sd_volume_select_siv_15_end[m] = np.nanstd(array_volume)
    area_select_siv_15_change[m] = 100. * (area_select_siv_15_end[m] - area_select_siv_15_start[m]) / area_select_siv_15_start[m]
    volume_select_siv_15_change[m] = 100. * (volume_select_siv_15_end[m] - volume_select_siv_15_start[m]) / volume_select_siv_15_start[m]

# Compute multi-model mean - good trend in SIA (15 best models) - Start
area_select_tsia_15_start = np.zeros(nmy)
volume_select_tsia_15_start = np.zeros(nmy)
sd_area_select_tsia_15_start = np.zeros(nmy)
sd_volume_select_tsia_15_start = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_hadgem3mm_start[m],area_cesm2waccm_start[m],area_fgoalsf3l_start[m],area_mpiesmhr_start[m],area_cesm2_start[m],area_mriesm_start[m],area_accessesm_start[m],area_bcccsm2mr_start[m],area_fioesm_start[m],area_gfdlesm4_start[m],area_ipslcm6alr_start[m],area_mpiesmlr_start[m],area_cnrmcm6hr_start[m],area_miroc6_start[m],area_gfdlcm4_start[m]]
        array_volume = [volume_hadgem3mm_start[m],volume_cesm2waccm_start[m],volume_fgoalsf3l_start[m],volume_mpiesmhr_start[m],volume_cesm2_start[m],volume_mriesm_start[m],volume_accessesm_start[m],volume_bcccsm2mr_start[m],volume_fioesm_start[m],volume_gfdlesm4_start[m],volume_ipslcm6alr_start[m],volume_mpiesmlr_start[m],volume_cnrmcm6hr_start[m],volume_miroc6_start[m],volume_gfdlcm4_start[m]]
    else:
        array_area = [area_hadgem3mm_start[m],area_cesm2waccm_start[m],area_fgoalsf3l_start[m],area_mpiesmhr_start[m],area_cesm2_start[m],area_mriesm_start[m],area_accessesm_start[m],area_bcccsm2mr_start[m],area_fioesm_start[m],area_gfdlesm4_start[m],area_ipslcm6alr_start[m],area_mpiesmlr_start[m],area_cnrmcm6hr_start[m],area_miroc6_start[m],area_hadgem3ll_start[m]]
        array_volume = [volume_hadgem3mm_start[m],volume_cesm2waccm_start[m],volume_fgoalsf3l_start[m],volume_mpiesmhr_start[m],volume_cesm2_start[m],volume_mriesm_start[m],volume_accessesm_start[m],volume_bcccsm2mr_start[m],volume_fioesm_start[m],volume_gfdlesm4_start[m],volume_ipslcm6alr_start[m],volume_mpiesmlr_start[m],volume_cnrmcm6hr_start[m],volume_miroc6_start[m],volume_hadgem3ll_start[m]]
    area_select_tsia_15_start[m] = np.nanmean(array_area)
    volume_select_tsia_15_start[m] = np.nanmean(array_volume)
    sd_area_select_tsia_15_start[m] = np.nanstd(array_area)
    sd_volume_select_tsia_15_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good trend in SIA (15 best models) - End
area_select_tsia_15_end = np.zeros(nmy)
volume_select_tsia_15_end = np.zeros(nmy)
sd_area_select_tsia_15_end = np.zeros(nmy)
sd_volume_select_tsia_15_end = np.zeros(nmy)
area_select_tsia_15_change = np.zeros(nmy)
volume_select_tsia_15_change = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_hadgem3mm_end[m],area_cesm2waccm_end[m],area_fgoalsf3l_end[m],area_mpiesmhr_end[m],area_cesm2_end[m],area_mriesm_end[m],area_accessesm_end[m],area_bcccsm2mr_end[m],area_fioesm_end[m],area_gfdlesm4_end[m],area_ipslcm6alr_end[m],area_mpiesmlr_end[m],area_cnrmcm6hr_end[m],area_miroc6_end[m],area_gfdlcm4_end[m]]
        array_volume = [volume_hadgem3mm_end[m],volume_cesm2waccm_end[m],volume_fgoalsf3l_end[m],volume_mpiesmhr_end[m],volume_cesm2_end[m],volume_mriesm_end[m],volume_accessesm_end[m],volume_bcccsm2mr_end[m],volume_fioesm_end[m],volume_gfdlesm4_end[m],volume_ipslcm6alr_end[m],volume_mpiesmlr_end[m],volume_cnrmcm6hr_end[m],volume_miroc6_end[m],volume_gfdlcm4_end[m]]
    else:
        array_area = [area_hadgem3mm_end[m],area_cesm2waccm_end[m],area_fgoalsf3l_end[m],area_mpiesmhr_end[m],area_cesm2_end[m],area_mriesm_end[m],area_accessesm_end[m],area_bcccsm2mr_end[m],area_fioesm_end[m],area_gfdlesm4_end[m],area_ipslcm6alr_end[m],area_mpiesmlr_end[m],area_cnrmcm6hr_end[m],area_miroc6_end[m],area_hadgem3ll_end[m]]
        array_volume = [volume_hadgem3mm_end[m],volume_cesm2waccm_end[m],volume_fgoalsf3l_end[m],volume_mpiesmhr_end[m],volume_cesm2_end[m],volume_mriesm_end[m],volume_accessesm_end[m],volume_bcccsm2mr_end[m],volume_fioesm_end[m],volume_gfdlesm4_end[m],volume_ipslcm6alr_end[m],volume_mpiesmlr_end[m],volume_cnrmcm6hr_end[m],volume_miroc6_end[m],volume_hadgem3ll_end[m]]
    area_select_tsia_15_end[m] = np.nanmean(array_area)
    volume_select_tsia_15_end[m] = np.nanmean(array_volume)
    sd_area_select_tsia_15_end[m] = np.nanstd(array_area)
    sd_volume_select_tsia_15_end[m] = np.nanstd(array_volume)
    area_select_tsia_15_change[m] = 100. * (area_select_tsia_15_end[m] - area_select_tsia_15_start[m]) / area_select_tsia_15_start[m]
    volume_select_tsia_15_change[m] = 100. * (volume_select_tsia_15_end[m] - volume_select_tsia_15_start[m]) / volume_select_tsia_15_start[m]

# Compute multi-model mean - good trend in SIV (15 best models) - Start
area_select_tsiv_15_start = np.zeros(nmy)
volume_select_tsiv_15_start = np.zeros(nmy)
sd_area_select_tsiv_15_start = np.zeros(nmy)
sd_volume_select_tsiv_15_start = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_ipslcm6alr_start[m],area_cesm2_start[m],area_fioesm_start[m],area_nesm3_start[m],area_noresm2mm_start[m],area_noresm2lm_start[m],area_mpiesmhr_start[m],area_miroc6_start[m],area_accesscm_start[m],area_gfdlcm4_start[m],area_mpiesmlr_start[m],area_fgoalsf3l_start[m],area_gfdlesm4_start[m],area_kiostesm_start[m],area_ecearth3_start[m]]
        array_volume = [volume_ipslcm6alr_start[m],volume_cesm2_start[m],volume_fioesm_start[m],volume_nesm3_start[m],volume_noresm2mm_start[m],volume_noresm2lm_start[m],volume_mpiesmhr_start[m],volume_miroc6_start[m],volume_accesscm_start[m],volume_gfdlcm4_start[m],volume_mpiesmlr_start[m],volume_fgoalsf3l_start[m],volume_gfdlesm4_start[m],volume_kiostesm_start[m],volume_ecearth3_start[m]]
    else:
        array_area = [area_ipslcm6alr_start[m],area_cesm2_start[m],area_fioesm_start[m],area_nesm3_start[m],area_noresm2mm_start[m],area_noresm2lm_start[m],area_mpiesmhr_start[m],area_miroc6_start[m],area_accesscm_start[m],area_mpiesmlr_start[m],area_fgoalsf3l_start[m],area_gfdlesm4_start[m],area_kiostesm_start[m],area_ecearth3_start[m],area_mriesm_start[m]]
        array_volume = [volume_ipslcm6alr_start[m],volume_cesm2_start[m],volume_fioesm_start[m],volume_nesm3_start[m],volume_noresm2mm_start[m],volume_noresm2lm_start[m],volume_mpiesmhr_start[m],volume_miroc6_start[m],volume_accesscm_start[m],volume_mpiesmlr_start[m],volume_fgoalsf3l_start[m],volume_gfdlesm4_start[m],volume_kiostesm_start[m],volume_ecearth3_start[m],volume_mriesm_start[m]]
    area_select_tsiv_15_start[m] = np.nanmean(array_area)
    volume_select_tsiv_15_start[m] = np.nanmean(array_volume)
    sd_area_select_tsiv_15_start[m] = np.nanstd(array_area)
    sd_volume_select_tsiv_15_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good trend in SIV (15 best models) - End
area_select_tsiv_15_end = np.zeros(nmy)
volume_select_tsiv_15_end = np.zeros(nmy)
sd_area_select_tsiv_15_end = np.zeros(nmy)
sd_volume_select_tsiv_15_end = np.zeros(nmy)
area_select_tsiv_15_change = np.zeros(nmy)
volume_select_tsiv_15_change = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_ipslcm6alr_end[m],area_cesm2_end[m],area_fioesm_end[m],area_nesm3_end[m],area_noresm2mm_end[m],area_noresm2lm_end[m],area_mpiesmhr_end[m],area_miroc6_end[m],area_accesscm_end[m],area_gfdlcm4_end[m],area_mpiesmlr_end[m],area_fgoalsf3l_end[m],area_gfdlesm4_end[m],area_kiostesm_end[m],area_ecearth3_end[m]]
        array_volume = [volume_ipslcm6alr_end[m],volume_cesm2_end[m],volume_fioesm_end[m],volume_nesm3_end[m],volume_noresm2mm_end[m],volume_noresm2lm_end[m],volume_mpiesmhr_end[m],volume_miroc6_end[m],volume_accesscm_end[m],volume_gfdlcm4_end[m],volume_mpiesmlr_end[m],volume_fgoalsf3l_end[m],volume_gfdlesm4_end[m],volume_kiostesm_end[m],volume_ecearth3_end[m]]
    else:
        array_area = [area_ipslcm6alr_end[m],area_cesm2_end[m],area_fioesm_end[m],area_nesm3_end[m],area_noresm2mm_end[m],area_noresm2lm_end[m],area_mpiesmhr_end[m],area_miroc6_end[m],area_accesscm_end[m],area_mpiesmlr_end[m],area_fgoalsf3l_end[m],area_gfdlesm4_end[m],area_kiostesm_end[m],area_ecearth3_end[m],area_mriesm_end[m]]
        array_volume = [volume_ipslcm6alr_end[m],volume_cesm2_end[m],volume_fioesm_end[m],volume_nesm3_end[m],volume_noresm2mm_end[m],volume_noresm2lm_end[m],volume_mpiesmhr_end[m],volume_miroc6_end[m],volume_accesscm_end[m],volume_mpiesmlr_end[m],volume_fgoalsf3l_end[m],volume_gfdlesm4_end[m],volume_kiostesm_end[m],volume_ecearth3_end[m],volume_mriesm_end[m]]
    area_select_tsiv_15_end[m] = np.nanmean(array_area)
    volume_select_tsiv_15_end[m] = np.nanmean(array_volume)
    sd_area_select_tsiv_15_end[m] = np.nanstd(array_area)
    sd_volume_select_tsiv_15_end[m] = np.nanstd(array_volume)
    area_select_tsiv_15_change[m] = 100. * (area_select_tsiv_15_end[m] - area_select_tsiv_15_start[m]) / area_select_tsiv_15_start[m]
    volume_select_tsiv_15_change[m] = 100. * (volume_select_tsiv_15_end[m] - volume_select_tsiv_15_start[m]) / volume_select_tsiv_15_start[m]

# Compute multi-model mean - good SIA variability (15 best models) - Start
area_select_varsia_15_start = np.zeros(nmy)
volume_select_varsia_15_start = np.zeros(nmy)
sd_area_select_varsia_15_start = np.zeros(nmy)
sd_volume_select_varsia_15_start = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_gfdlcm4_start[m],area_accessesm_start[m],area_inmcm50_start[m],area_noresm2mm_start[m],area_gfdlesm4_start[m],area_canesm5canoe_start[m],area_cnrmcm6hr_start[m],area_accesscm_start[m],area_hadgem3ll_start[m],area_nesm3_start[m],area_fgoalsf3l_start[m],area_fgoalsg3_start[m],area_inmcm48_start[m],area_bcccsm2mr_start[m],area_ecearth3_start[m]]
        array_volume = [volume_gfdlcm4_start[m],volume_accessesm_start[m],volume_noresm2mm_start[m],volume_gfdlesm4_start[m],volume_cnrmcm6hr_start[m],volume_accesscm_start[m],volume_hadgem3ll_start[m],volume_nesm3_start[m],volume_fgoalsf3l_start[m],volume_bcccsm2mr_start[m],volume_ecearth3_start[m],volume_ecearth3veg_start[m],volume_cesm2waccm_start[m],volume_kiostesm_start[m],volume_camscsm_start[m]]
    else:
        array_area = [area_accessesm_start[m],area_inmcm50_start[m],area_noresm2mm_start[m],area_gfdlesm4_start[m],area_canesm5canoe_start[m],area_cnrmcm6hr_start[m],area_accesscm_start[m],area_hadgem3ll_start[m],area_nesm3_start[m],area_fgoalsf3l_start[m],area_fgoalsg3_start[m],area_inmcm48_start[m],area_bcccsm2mr_start[m],area_ecearth3_start[m],area_mriesm_start[m]]
        array_volume = [volume_accessesm_start[m],volume_noresm2mm_start[m],volume_gfdlesm4_start[m],volume_cnrmcm6hr_start[m],volume_accesscm_start[m],volume_hadgem3ll_start[m],volume_nesm3_start[m],volume_fgoalsf3l_start[m],volume_bcccsm2mr_start[m],volume_ecearth3_start[m],volume_ecearth3veg_start[m],volume_cesm2waccm_start[m],volume_kiostesm_start[m],volume_camscsm_start[m],volume_mriesm_start[m]] 
    area_select_varsia_15_start[m] = np.nanmean(array_area)
    volume_select_varsia_15_start[m] = np.nanmean(array_volume)
    sd_area_select_varsia_15_start[m] = np.nanstd(array_area)
    sd_volume_select_varsia_15_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good SIA variability (15 best models) - End
area_select_varsia_15_end = np.zeros(nmy)
volume_select_varsia_15_end = np.zeros(nmy)
sd_area_select_varsia_15_end = np.zeros(nmy)
sd_volume_select_varsia_15_end = np.zeros(nmy)
area_select_varsia_15_change = np.zeros(nmy)
volume_select_varsia_15_change = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_gfdlcm4_end[m],area_accessesm_end[m],area_inmcm50_end[m],area_noresm2mm_end[m],area_gfdlesm4_end[m],area_canesm5canoe_end[m],area_cnrmcm6hr_end[m],area_accesscm_end[m],area_hadgem3ll_end[m],area_nesm3_end[m],area_fgoalsf3l_end[m],area_fgoalsg3_end[m],area_inmcm48_end[m],area_bcccsm2mr_end[m],area_ecearth3_end[m]]
        array_volume = [volume_gfdlcm4_end[m],volume_accessesm_end[m],volume_noresm2mm_end[m],volume_gfdlesm4_end[m],volume_cnrmcm6hr_end[m],volume_accesscm_end[m],volume_hadgem3ll_end[m],volume_nesm3_end[m],volume_fgoalsf3l_end[m],volume_bcccsm2mr_end[m],volume_ecearth3_end[m],volume_ecearth3veg_end[m],volume_cesm2waccm_end[m],volume_kiostesm_end[m],volume_camscsm_end[m]]
    else:
        array_area = [area_accessesm_end[m],area_inmcm50_end[m],area_noresm2mm_end[m],area_gfdlesm4_end[m],area_canesm5canoe_end[m],area_cnrmcm6hr_end[m],area_accesscm_end[m],area_hadgem3ll_end[m],area_nesm3_end[m],area_fgoalsf3l_end[m],area_fgoalsg3_end[m],area_inmcm48_end[m],area_bcccsm2mr_end[m],area_ecearth3_end[m],area_mriesm_end[m]]
        array_volume = [volume_accessesm_end[m],volume_noresm2mm_end[m],volume_gfdlesm4_end[m],volume_cnrmcm6hr_end[m],volume_accesscm_end[m],volume_hadgem3ll_end[m],volume_nesm3_end[m],volume_fgoalsf3l_end[m],volume_bcccsm2mr_end[m],volume_ecearth3_end[m],volume_ecearth3veg_end[m],volume_cesm2waccm_end[m],volume_kiostesm_end[m],volume_camscsm_end[m],volume_mriesm_end[m]] 
    area_select_varsia_15_end[m] = np.nanmean(array_area)
    volume_select_varsia_15_end[m] = np.nanmean(array_volume)
    sd_area_select_varsia_15_end[m] = np.nanstd(array_area)
    sd_volume_select_varsia_15_end[m] = np.nanstd(array_volume)
    area_select_varsia_15_change[m] = 100. * (area_select_varsia_15_end[m] - area_select_varsia_15_start[m]) / area_select_varsia_15_start[m]
    volume_select_varsia_15_change[m] = 100. * (volume_select_varsia_15_end[m] - volume_select_varsia_15_start[m]) / volume_select_varsia_15_start[m]

# Compute multi-model mean - good SIV variability (15 best models) - Start
area_select_varsiv_15_start = np.zeros(nmy)
volume_select_varsiv_15_start = np.zeros(nmy)
sd_area_select_varsiv_15_start = np.zeros(nmy)
sd_volume_select_varsiv_15_start = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_gfdlcm4_start[m],area_canesm5_start[m],area_noresm2lm_start[m],area_noresm2mm_start[m],area_fioesm_start[m],area_hadgem3mm_start[m],area_cesm2waccm_start[m],area_miroces2l_start[m],area_accesscm_start[m],area_gfdlesm4_start[m],area_cesm2_start[m],area_ecearth3veg_start[m],area_ecearth3_start[m],area_fgoalsf3l_start[m],area_camscsm_start[m]]
        array_volume = [volume_gfdlcm4_start[m],volume_canesm5_start[m],volume_noresm2lm_start[m],volume_noresm2mm_start[m],volume_fioesm_start[m],volume_hadgem3mm_start[m],volume_cesm2waccm_start[m],volume_miroces2l_start[m],volume_accesscm_start[m],volume_gfdlesm4_start[m],volume_cesm2_start[m],volume_ecearth3veg_start[m],volume_ecearth3_start[m],volume_fgoalsf3l_start[m],volume_camscsm_start[m]]
    else:
        array_area = [area_canesm5_start[m],area_noresm2lm_start[m],area_noresm2mm_start[m],area_fioesm_start[m],area_hadgem3mm_start[m],area_cesm2waccm_start[m],area_miroces2l_start[m],area_accesscm_start[m],area_gfdlesm4_start[m],area_cesm2_start[m],area_ecearth3veg_start[m],area_ecearth3_start[m],area_fgoalsf3l_start[m],area_camscsm_start[m],area_ukesmll_start[m]]
        array_volume = [volume_canesm5_start[m],volume_noresm2lm_start[m],volume_noresm2mm_start[m],volume_fioesm_start[m],volume_hadgem3mm_start[m],volume_cesm2waccm_start[m],volume_miroces2l_start[m],volume_accesscm_start[m],volume_gfdlesm4_start[m],volume_cesm2_start[m],volume_ecearth3veg_start[m],volume_ecearth3_start[m],volume_fgoalsf3l_start[m],volume_camscsm_start[m],volume_ukesmll_start[m]]
    area_select_varsiv_15_start[m] = np.nanmean(array_area)
    volume_select_varsiv_15_start[m] = np.nanmean(array_volume)
    sd_area_select_varsiv_15_start[m] = np.nanstd(array_area)
    sd_volume_select_varsiv_15_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good SIV variability (15 best models) - End
area_select_varsiv_15_end = np.zeros(nmy)
volume_select_varsiv_15_end = np.zeros(nmy)
sd_area_select_varsiv_15_end = np.zeros(nmy)
sd_volume_select_varsiv_15_end = np.zeros(nmy)
area_select_varsiv_15_change = np.zeros(nmy)
volume_select_varsiv_15_change = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_gfdlcm4_end[m],area_canesm5_end[m],area_noresm2lm_end[m],area_noresm2mm_end[m],area_fioesm_end[m],area_hadgem3mm_end[m],area_cesm2waccm_end[m],area_miroces2l_end[m],area_accesscm_end[m],area_gfdlesm4_end[m],area_cesm2_end[m],area_ecearth3veg_end[m],area_ecearth3_end[m],area_fgoalsf3l_end[m],area_camscsm_end[m]]
        array_volume = [volume_gfdlcm4_end[m],volume_canesm5_end[m],volume_noresm2lm_end[m],volume_noresm2mm_end[m],volume_fioesm_end[m],volume_hadgem3mm_end[m],volume_cesm2waccm_end[m],volume_miroces2l_end[m],volume_accesscm_end[m],volume_gfdlesm4_end[m],volume_cesm2_end[m],volume_ecearth3veg_end[m],volume_ecearth3_end[m],volume_fgoalsf3l_end[m],volume_camscsm_end[m]]
    else:
        array_area = [area_canesm5_end[m],area_noresm2lm_end[m],area_noresm2mm_end[m],area_fioesm_end[m],area_hadgem3mm_end[m],area_cesm2waccm_end[m],area_miroces2l_end[m],area_accesscm_end[m],area_gfdlesm4_end[m],area_cesm2_end[m],area_ecearth3veg_end[m],area_ecearth3_end[m],area_fgoalsf3l_end[m],area_camscsm_end[m],area_ukesmll_end[m]]
        array_volume = [volume_canesm5_end[m],volume_noresm2lm_end[m],volume_noresm2mm_end[m],volume_fioesm_end[m],volume_hadgem3mm_end[m],volume_cesm2waccm_end[m],volume_miroces2l_end[m],volume_accesscm_end[m],volume_gfdlesm4_end[m],volume_cesm2_end[m],volume_ecearth3veg_end[m],volume_ecearth3_end[m],volume_fgoalsf3l_end[m],volume_camscsm_end[m],volume_ukesmll_end[m]]
    area_select_varsiv_15_end[m] = np.nanmean(array_area)
    volume_select_varsiv_15_end[m] = np.nanmean(array_volume)
    sd_area_select_varsiv_15_end[m] = np.nanstd(array_area)
    sd_volume_select_varsiv_15_end[m] = np.nanstd(array_volume)
    area_select_varsiv_15_change[m] = 100. * (area_select_varsiv_15_end[m] - area_select_varsiv_15_start[m]) / area_select_varsiv_15_start[m]
    volume_select_varsiv_15_change[m] = 100. * (volume_select_varsiv_15_end[m] - volume_select_varsiv_15_start[m]) / volume_select_varsiv_15_start[m]

# Compute multi-model mean - good AOHT 26N and 57N (best 8 models) - Start
area_select_aoht_8_start = np.zeros(nmy)
volume_select_aoht_8_start = np.zeros(nmy)
sd_area_select_aoht_8_start = np.zeros(nmy)
sd_volume_select_aoht_8_start = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_hadgem3mm_start[m],area_ecearth3veg_start[m],area_mriesm_start[m],area_mpiesmlr_start[m],area_ecearth3_start[m],area_hadgem3ll_start[m],area_ukesmll_start[m],area_cnrmcm6_start[m]]
    array_volume = [volume_hadgem3mm_start[m],volume_ecearth3veg_start[m],volume_mriesm_start[m],volume_mpiesmlr_start[m],volume_ecearth3_start[m],volume_hadgem3ll_start[m],volume_ukesmll_start[m],volume_cnrmcm6_start[m]]
    area_select_aoht_8_start[m] = np.nanmean(array_area)
    volume_select_aoht_8_start[m] = np.nanmean(array_volume)
    sd_area_select_aoht_8_start[m] = np.nanstd(array_area)
    sd_volume_select_aoht_8_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good AOHT 26N and 57N (best 8 models) - End
area_select_aoht_8_end = np.zeros(nmy)
volume_select_aoht_8_end = np.zeros(nmy)
sd_area_select_aoht_8_end = np.zeros(nmy)
sd_volume_select_aoht_8_end = np.zeros(nmy)
area_select_aoht_8_change = np.zeros(nmy)
volume_select_aoht_8_change = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_hadgem3mm_end[m],area_ecearth3veg_end[m],area_mriesm_end[m],area_mpiesmlr_end[m],area_ecearth3_end[m],area_hadgem3ll_end[m],area_ukesmll_end[m],area_cnrmcm6_end[m]]
    array_volume = [volume_hadgem3mm_end[m],volume_ecearth3veg_end[m],volume_mriesm_end[m],volume_mpiesmlr_end[m],volume_ecearth3_end[m],volume_hadgem3ll_end[m],volume_ukesmll_end[m],volume_cnrmcm6_end[m]]
    area_select_aoht_8_end[m] = np.nanmean(array_area)
    volume_select_aoht_8_end[m] = np.nanmean(array_volume)
    sd_area_select_aoht_8_end[m] = np.nanstd(array_area)
    sd_volume_select_aoht_8_end[m] = np.nanstd(array_volume)
    area_select_aoht_8_change[m] = 100. * (area_select_aoht_8_end[m] - area_select_aoht_8_start[m]) / area_select_aoht_8_start[m]
    volume_select_aoht_8_change[m] = 100. * (volume_select_aoht_8_end[m] - volume_select_aoht_8_start[m]) / volume_select_aoht_8_start[m]

# Compute multi-model mean - good AOHT 70N and POHT 60N (best 8 models) - Start
area_select_apoht_8_start = np.zeros(nmy)
volume_select_apoht_8_start = np.zeros(nmy)
sd_area_select_apoht_8_start = np.zeros(nmy)
sd_volume_select_apoht_8_start = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_mpiesmlr_start[m],area_ukesmll_start[m],area_mpiesmhr_start[m],area_hadgem3ll_start[m],area_hadgem3mm_start[m],area_ipslcm6alr_start[m],area_mriesm_start[m],area_canesm5_start[m]]
    array_volume = [volume_mpiesmlr_start[m],volume_ukesmll_start[m],volume_mpiesmhr_start[m],volume_hadgem3ll_start[m],volume_hadgem3mm_start[m],volume_ipslcm6alr_start[m],volume_mriesm_start[m],volume_canesm5_start[m]]
    area_select_apoht_8_start[m] = np.nanmean(array_area)
    volume_select_apoht_8_start[m] = np.nanmean(array_volume)
    sd_area_select_apoht_8_start[m] = np.nanstd(array_area)
    sd_volume_select_apoht_8_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good AOHT 70N and POHT 60N (best 8 models) - End
area_select_apoht_8_end = np.zeros(nmy)
volume_select_apoht_8_end = np.zeros(nmy)
sd_area_select_apoht_8_end = np.zeros(nmy)
sd_volume_select_apoht_8_end = np.zeros(nmy)
area_select_apoht_8_change = np.zeros(nmy)
volume_select_apoht_8_change = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_mpiesmlr_end[m],area_ukesmll_end[m],area_mpiesmhr_end[m],area_hadgem3ll_end[m],area_hadgem3mm_end[m],area_ipslcm6alr_end[m],area_mriesm_end[m],area_canesm5_end[m]]
    array_volume = [volume_mpiesmlr_end[m],volume_ukesmll_end[m],volume_mpiesmhr_end[m],volume_hadgem3ll_end[m],volume_hadgem3mm_end[m],volume_ipslcm6alr_end[m],volume_mriesm_end[m],volume_canesm5_end[m]]
    area_select_apoht_8_end[m] = np.nanmean(array_area)
    volume_select_apoht_8_end[m] = np.nanmean(array_volume)
    sd_area_select_apoht_8_end[m] = np.nanstd(array_area)
    sd_volume_select_apoht_8_end[m] = np.nanstd(array_volume)
    area_select_apoht_8_change[m] = 100. * (area_select_apoht_8_end[m] - area_select_apoht_8_start[m]) / area_select_apoht_8_start[m]
    volume_select_apoht_8_change[m] = 100. * (volume_select_apoht_8_end[m] - volume_select_apoht_8_start[m]) / volume_select_apoht_8_start[m]

# Compute multi-model mean - good AOHT 26N+57N and SIA (among the best 15 SIA models and best 8 OHT models) - Start
area_select_ohtsia1_start = np.zeros(nmy)
volume_select_ohtsia1_start = np.zeros(nmy)
sd_area_select_ohtsia1_start = np.zeros(nmy)
sd_volume_select_ohtsia1_start = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_hadgem3mm_start[m],area_ecearth3veg_start[m],area_mpiesmlr_start[m],area_ecearth3_start[m],area_hadgem3ll_start[m],area_cnrmcm6_start[m]]
    array_volume = [volume_hadgem3mm_start[m],volume_ecearth3veg_start[m],volume_mpiesmlr_start[m],volume_ecearth3_start[m],volume_hadgem3ll_start[m],volume_cnrmcm6_start[m]]
    area_select_ohtsia1_start[m] = np.nanmean(array_area)
    volume_select_ohtsia1_start[m] = np.nanmean(array_volume)
    sd_area_select_ohtsia1_start[m] = np.nanstd(array_area)
    sd_volume_select_ohtsia1_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good AOHT 26N+57N and SIA (among the best 15 SIA models and best 8 OHT models) - End
area_select_ohtsia1_end = np.zeros(nmy)
volume_select_ohtsia1_end = np.zeros(nmy)
sd_area_select_ohtsia1_end = np.zeros(nmy)
sd_volume_select_ohtsia1_end = np.zeros(nmy)
area_select_ohtsia1_change = np.zeros(nmy)
volume_select_ohtsia1_change = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_hadgem3mm_end[m],area_ecearth3veg_end[m],area_mpiesmlr_end[m],area_ecearth3_end[m],area_hadgem3ll_end[m],area_cnrmcm6_end[m]]
    array_volume = [volume_hadgem3mm_end[m],volume_ecearth3veg_end[m],volume_mpiesmlr_end[m],volume_ecearth3_end[m],volume_hadgem3ll_end[m],volume_cnrmcm6_end[m]]
    area_select_ohtsia1_end[m] = np.nanmean(array_area)
    volume_select_ohtsia1_end[m] = np.nanmean(array_volume)
    sd_area_select_ohtsia1_end[m] = np.nanstd(array_area)
    sd_volume_select_ohtsia1_end[m] = np.nanstd(array_volume)
    area_select_ohtsia1_change[m] = 100. * (area_select_ohtsia1_end[m] - area_select_ohtsia1_start[m]) / area_select_ohtsia1_start[m]
    volume_select_ohtsia1_change[m] = 100. * (volume_select_ohtsia1_end[m] - volume_select_ohtsia1_start[m]) / volume_select_ohtsia1_start[m]

# Compute multi-model mean - good AOHT 70N + POHT 60N and SIA (among the best 15 SIA models and best 8 OHT models) - Start
area_select_ohtsia2_start = np.zeros(nmy)
volume_select_ohtsia2_start = np.zeros(nmy)
sd_area_select_ohtsia2_start = np.zeros(nmy)
sd_volume_select_ohtsia2_start = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_mpiesmlr_start[m],area_hadgem3ll_start[m],area_hadgem3mm_start[m],area_ipslcm6alr_start[m],area_canesm5_start[m]]
    array_volume = [volume_mpiesmlr_start[m],volume_hadgem3ll_start[m],volume_hadgem3mm_start[m],volume_ipslcm6alr_start[m],volume_canesm5_start[m]]
    area_select_ohtsia2_start[m] = np.nanmean(array_area)
    volume_select_ohtsia2_start[m] = np.nanmean(array_volume)
    sd_area_select_ohtsia2_start[m] = np.nanstd(array_area)
    sd_volume_select_ohtsia2_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good AOHT 70N + POHT 60N and SIA (among the best 15 SIA models and best 8 OHT models) - End
area_select_ohtsia2_end = np.zeros(nmy)
volume_select_ohtsia2_end = np.zeros(nmy)
sd_area_select_ohtsia2_end = np.zeros(nmy)
sd_volume_select_ohtsia2_end = np.zeros(nmy)
area_select_ohtsia2_change = np.zeros(nmy)
volume_select_ohtsia2_change = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_mpiesmlr_end[m],area_hadgem3ll_end[m],area_hadgem3mm_end[m],area_ipslcm6alr_end[m],area_canesm5_end[m]]
    array_volume = [volume_mpiesmlr_end[m],volume_hadgem3ll_end[m],volume_hadgem3mm_end[m],volume_ipslcm6alr_end[m],volume_canesm5_end[m]]
    area_select_ohtsia2_end[m] = np.nanmean(array_area)
    volume_select_ohtsia2_end[m] = np.nanmean(array_volume)
    sd_area_select_ohtsia2_end[m] = np.nanstd(array_area)
    sd_volume_select_ohtsia2_end[m] = np.nanstd(array_volume)
    area_select_ohtsia2_change[m] = 100. * (area_select_ohtsia2_end[m] - area_select_ohtsia2_start[m]) / area_select_ohtsia2_start[m]
    volume_select_ohtsia2_change[m] = 100. * (volume_select_ohtsia2_end[m] - volume_select_ohtsia2_start[m]) / volume_select_ohtsia2_start[m]

# Compute multi-model mean - good AOHT 26N+57N and SIV (among the best 15 SIV models and best 8 OHT models) - Start
area_select_ohtsiv1_start = np.zeros(nmy)
volume_select_ohtsiv1_start = np.zeros(nmy)
sd_area_select_ohtsiv1_start = np.zeros(nmy)
sd_volume_select_ohtsiv1_start = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_hadgem3mm_start[m],area_mriesm_start[m],area_mpiesmlr_start[m]]
    array_volume = [volume_hadgem3mm_start[m],volume_mriesm_start[m],volume_mpiesmlr_start[m]]
    area_select_ohtsiv1_start[m] = np.nanmean(array_area)
    volume_select_ohtsiv1_start[m] = np.nanmean(array_volume)
    sd_area_select_ohtsiv1_start[m] = np.nanstd(array_area)
    sd_volume_select_ohtsiv1_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good AOHT 26N+57N and SIV (among the best 15 SIV models and best 8 OHT models) - End
area_select_ohtsiv1_end = np.zeros(nmy)
volume_select_ohtsiv1_end = np.zeros(nmy)
sd_area_select_ohtsiv1_end = np.zeros(nmy)
sd_volume_select_ohtsiv1_end = np.zeros(nmy)
area_select_ohtsiv1_change = np.zeros(nmy)
volume_select_ohtsiv1_change = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_hadgem3mm_end[m],area_mriesm_end[m],area_mpiesmlr_end[m]]
    array_volume = [volume_hadgem3mm_end[m],volume_mriesm_end[m],volume_mpiesmlr_end[m]]
    area_select_ohtsiv1_end[m] = np.nanmean(array_area)
    volume_select_ohtsiv1_end[m] = np.nanmean(array_volume)
    sd_area_select_ohtsiv1_end[m] = np.nanstd(array_area)
    sd_volume_select_ohtsiv1_end[m] = np.nanstd(array_volume)
    area_select_ohtsiv1_change[m] = 100. * (area_select_ohtsiv1_end[m] - area_select_ohtsiv1_start[m]) / area_select_ohtsiv1_start[m]
    volume_select_ohtsiv1_change[m] = 100. * (volume_select_ohtsiv1_end[m] - volume_select_ohtsiv1_start[m]) / volume_select_ohtsiv1_start[m]

# Compute multi-model mean - good AOHT 70N + POHT 60N and SIV (among the best 15 SIV models and best 8 OHT models) - Start
area_select_ohtsiv2_start = np.zeros(nmy)
volume_select_ohtsiv2_start = np.zeros(nmy)
sd_area_select_ohtsiv2_start = np.zeros(nmy)
sd_volume_select_ohtsiv2_start = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_mpiesmlr_start[m],area_mpiesmhr_start[m],area_hadgem3mm_start[m],area_ipslcm6alr_start[m],area_mriesm_start[m]]
    array_volume = [volume_mpiesmlr_start[m],volume_mpiesmhr_start[m],volume_hadgem3mm_start[m],volume_ipslcm6alr_start[m],volume_mriesm_start[m]]
    area_select_ohtsiv2_start[m] = np.nanmean(array_area)
    volume_select_ohtsiv2_start[m] = np.nanmean(array_volume)
    sd_area_select_ohtsiv2_start[m] = np.nanstd(array_area)
    sd_volume_select_ohtsiv2_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - good AOHT 70N + POHT 60N and SIV (among the best 15 SIV models and best 8 OHT models) - End
area_select_ohtsiv2_end = np.zeros(nmy)
volume_select_ohtsiv2_end = np.zeros(nmy)
sd_area_select_ohtsiv2_end = np.zeros(nmy)
sd_volume_select_ohtsiv2_end = np.zeros(nmy)
area_select_ohtsiv2_change = np.zeros(nmy)
volume_select_ohtsiv2_change = np.zeros(nmy)
for m in np.arange(nmy):
    array_area = [area_mpiesmlr_end[m],area_mpiesmhr_end[m],area_hadgem3mm_end[m],area_ipslcm6alr_end[m],area_mriesm_end[m]]
    array_volume = [volume_mpiesmlr_end[m],volume_mpiesmhr_end[m],volume_hadgem3mm_end[m],volume_ipslcm6alr_end[m],volume_mriesm_end[m]]
    area_select_ohtsiv2_end[m] = np.nanmean(array_area)
    volume_select_ohtsiv2_end[m] = np.nanmean(array_volume)
    sd_area_select_ohtsiv2_end[m] = np.nanstd(array_area)
    sd_volume_select_ohtsiv2_end[m] = np.nanstd(array_volume)
    area_select_ohtsiv2_change[m] = 100. * (area_select_ohtsiv2_end[m] - area_select_ohtsiv2_start[m]) / area_select_ohtsiv2_start[m]
    volume_select_ohtsiv2_change[m] = 100. * (volume_select_ohtsiv2_end[m] - volume_select_ohtsiv2_start[m]) / volume_select_ohtsiv2_start[m]

# Compute multi-model mean - >=5 members - Start
area_select_members_start = np.zeros(nmy)
volume_select_members_start = np.zeros(nmy)
sd_area_select_members_start = np.zeros(nmy)
sd_volume_select_members_start = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_canesm5_start[m],area_cnrmcm6_start[m],area_mpiesmlr_start[m],area_ecearth3_start[m],area_ecearth3veg_start[m],area_ipslcm6alr_start[m],area_miroc6_start[m],area_ukesmll_start[m],area_cesm2_start[m],area_cesm2waccm_start[m]]
        array_volume = [volume_canesm5_start[m],volume_cnrmcm6_start[m],volume_mpiesmlr_start[m],volume_ecearth3_start[m],volume_ecearth3veg_start[m],volume_ipslcm6alr_start[m],volume_miroc6_start[m],volume_ukesmll_start[m],volume_cesm2_start[m],volume_cesm2waccm_start[m]]
    elif experiment == 'ssp126':
        array_area = [area_canesm5_start[m],area_cnrmcm6_start[m],area_mpiesmlr_start[m],area_ecearth3_start[m],area_ecearth3veg_start[m],area_ipslcm6alr_start[m],area_miroc6_start[m],area_ukesmll_start[m],area_cesm2_start[m]]
        array_volume = [volume_canesm5_start[m],volume_cnrmcm6_start[m],volume_mpiesmlr_start[m],volume_ecearth3_start[m],volume_ecearth3veg_start[m],volume_ipslcm6alr_start[m],volume_miroc6_start[m],volume_ukesmll_start[m],volume_cesm2_start[m]]
    area_select_members_start[m] = np.nanmean(array_area)
    volume_select_members_start[m] = np.nanmean(array_volume)
    sd_area_select_members_start[m] = np.nanstd(array_area)
    sd_volume_select_members_start[m] = np.nanstd(array_volume)

# Compute multi-model mean - >=5 members - End
area_select_members_end = np.zeros(nmy)
volume_select_members_end = np.zeros(nmy)
sd_area_select_members_end = np.zeros(nmy)
sd_volume_select_members_end = np.zeros(nmy)
area_select_members_change = np.zeros(nmy)
volume_select_members_change = np.zeros(nmy)
for m in np.arange(nmy):
    if experiment == 'ssp585':
        array_area = [area_canesm5_end[m],area_cnrmcm6_end[m],area_mpiesmlr_end[m],area_ecearth3_end[m],area_ecearth3veg_end[m],area_ipslcm6alr_end[m],area_miroc6_end[m],area_ukesmll_end[m],area_cesm2_end[m],area_cesm2waccm_end[m]]
        array_volume = [volume_canesm5_end[m],volume_cnrmcm6_end[m],volume_mpiesmlr_end[m],volume_ecearth3_end[m],volume_ecearth3veg_end[m],volume_ipslcm6alr_end[m],volume_miroc6_end[m],volume_ukesmll_end[m],volume_cesm2_end[m],volume_cesm2waccm_end[m]]
    elif experiment == 'ssp126':
        array_area = [area_canesm5_end[m],area_cnrmcm6_end[m],area_mpiesmlr_end[m],area_ecearth3_end[m],area_ecearth3veg_end[m],area_ipslcm6alr_end[m],area_miroc6_end[m],area_ukesmll_end[m],area_cesm2_end[m]]
        array_volume = [volume_canesm5_end[m],volume_cnrmcm6_end[m],volume_mpiesmlr_end[m],volume_ecearth3_end[m],volume_ecearth3veg_end[m],volume_ipslcm6alr_end[m],volume_miroc6_end[m],volume_ukesmll_end[m],volume_cesm2_end[m]]
    area_select_members_end[m] = np.nanmean(array_area)
    volume_select_members_end[m] = np.nanmean(array_volume)
    sd_area_select_members_end[m] = np.nanstd(array_area)
    sd_volume_select_members_end[m] = np.nanstd(array_volume)
    area_select_members_change[m] = 100. * (area_select_members_end[m] - area_select_members_start[m]) / area_select_members_start[m]
    volume_select_members_change[m] = 100. * (volume_select_members_end[m] - volume_select_members_start[m]) / volume_select_members_start[m]

# Scatter plots
fig,ax = plt.subplots(2,2,figsize=(16,15))
fig.subplots_adjust(left=0.1,bottom=0.2,right=0.95,top=0.95,wspace=0.3,hspace=0.3)

# Change in March sea-ice area (2096-2100 vs 2015-2019)
month = 2
ax[0,0].plot(area_mmm_start[month],area_mmm_end[month],'o',color='black',markersize=16)
ax[0,0].errorbar(area_mmm_start[month],area_mmm_end[month],xerr=sd_area_mmm_start[month],yerr=sd_area_mmm_end[month],fmt='ko',capsize=5)
ax[0,0].plot(area_select_sia_15_start[month],area_select_sia_15_end[month],'o',color='blue',markersize=14)
ax[0,0].plot(area_select_siv_15_start[month],area_select_siv_15_end[month],'o',color='green',markersize=14)
ax[0,0].plot(area_select_varsia_15_start[month],area_select_varsia_15_end[month],'X',color='blue',markersize=14)
ax[0,0].plot(area_select_varsiv_15_start[month],area_select_varsiv_15_end[month],'X',color='green',markersize=14)
ax[0,0].plot(area_select_tsia_15_start[month],area_select_tsia_15_end[month],'P',color='blue',markersize=14)
ax[0,0].plot(area_select_tsiv_15_start[month],area_select_tsiv_15_end[month],'P',color='green',markersize=14)
ax[0,0].plot(area_select_aoht_8_start[month],area_select_aoht_8_end[month],'o',color='red',markersize=14)
ax[0,0].plot(area_select_apoht_8_start[month],area_select_apoht_8_end[month],'X',color='red',markersize=14)
ax[0,0].plot(area_select_ohtsia1_start[month],area_select_ohtsia1_end[month],'o',color='orange',markersize=14)
ax[0,0].plot(area_select_ohtsia2_start[month],area_select_ohtsia2_end[month],'X',color='orange',markersize=14)
ax[0,0].plot(area_select_ohtsiv1_start[month],area_select_ohtsiv1_end[month],'o',color='gray',markersize=14)
ax[0,0].plot(area_select_ohtsiv2_start[month],area_select_ohtsiv2_end[month],'X',color='gray',markersize=14)
ax[0,0].plot(area_select_members_start[month],area_select_members_end[month],'o',color='purple',markersize=14)
ax[0,0].plot(area_select_random_start[month],area_select_random_end[month],'o',color='cyan',markersize=10)
ax[0,0].errorbar(area_select_random_start[month],area_select_random_end[month],xerr=sd_area_select_random_start[month],yerr=sd_area_select_random_end[month],fmt='o',color='cyan',capsize=5)
ax[0,0].axvline(x=area_obs_start[month],color='black',linestyle='--',linewidth=2)
if experiment == 'ssp585':
    ax[0,0].annotate(str(int(np.round(area_mmm_change[month])))+'%',(area_mmm_start[month]+0.05,area_mmm_end[month]+0.2),color='black',fontsize=16)
    ax[0,0].annotate(str(int(np.round(area_select_sia_15_change[month])))+'%',(area_select_sia_15_start[month]-0.1,area_select_sia_15_end[month]-0.6),color='blue',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_siv_15_change[month])))+'%',(area_select_siv_15_start[month]+0.05,area_select_sia_15_end[month]-0.3),color='green',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_varsia_15_change[month])))+'%',(area_select_varsia_15_start[month]+0.05,area_select_varsia_15_end[month]+0.2),color='blue',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_varsiv_15_change[month])))+'%',(area_select_varsiv_15_start[month]+0.1,area_select_varsiv_15_end[month]-0.25),color='green',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_tsia_15_change[month])))+'%',(area_select_tsia_15_start[month]-0.3,area_select_tsia_15_end[month]+0.1),color='blue',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_tsiv_15_change[month])))+'%',(area_select_tsiv_15_start[month]-0.3,area_select_tsiv_15_end[month]+0.1),color='green',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_aoht_8_change[month])))+'%',(area_select_aoht_8_start[month]+0.08,area_select_aoht_8_end[month]-0.2),color='red',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_apoht_8_change[month])))+'%',(area_select_apoht_8_start[month]+0.05,area_select_apoht_8_end[month]+0.1),color='red',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_ohtsia1_change[month])))+'%',(area_select_ohtsia1_start[month]-0.1,area_select_ohtsia1_end[month]+0.25),color='orange',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_ohtsia2_change[month])))+'%',(area_select_ohtsia2_start[month]+0.05,area_select_ohtsia2_end[month]+0.2),color='orange',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_ohtsiv1_change[month])))+'%',(area_select_ohtsiv1_start[month]+0.05,area_select_ohtsiv1_end[month]+0.2),color='gray',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_ohtsiv2_change[month])))+'%',(area_select_ohtsiv2_start[month]-0.2,area_select_ohtsiv2_end[month]-0.5),color='gray',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_members_change[month])))+'%',(area_select_members_start[month]-0.3,area_select_members_end[month]-0.5),color='purple',fontsize=14)
else:
    ax[0,0].annotate(str(int(np.round(area_mmm_change[month])))+'%',(area_mmm_start[month]+0.05,area_mmm_end[month]-0.2),color='black',fontsize=16)
    ax[0,0].annotate(str(int(np.round(area_select_sia_15_change[month])))+'%',(area_select_sia_15_start[month]+0.07,area_select_sia_15_end[month]-0.05),color='blue',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_siv_15_change[month])))+'%',(area_select_siv_15_start[month]-0.05,area_select_siv_15_end[month]+0.1),color='green',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_varsia_15_change[month])))+'%',(area_select_varsia_15_start[month]+0.05,area_select_varsia_15_end[month]+0.1),color='blue',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_varsiv_15_change[month])))+'%',(area_select_varsiv_15_start[month]+0.05,area_select_varsiv_15_end[month]+0.1),color='green',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_tsia_15_change[month])))+'%',(area_select_tsia_15_start[month]-0.2,area_select_tsia_15_end[month]+0.1),color='blue',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_tsiv_15_change[month])))+'%',(area_select_tsiv_15_start[month]-0.1,area_select_tsiv_15_end[month]+0.1),color='green',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_aoht_8_change[month])))+'%',(area_select_aoht_8_start[month]-0.3,area_select_aoht_8_end[month]+0.05),color='red',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_apoht_8_change[month])))+'%',(area_select_apoht_8_start[month]+0.08,area_select_apoht_8_end[month]-0.05),color='red',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_ohtsia1_change[month])))+'%',(area_select_ohtsia1_start[month]+0.08,area_select_ohtsia1_end[month]),color='orange',fontsize=13)
    ax[0,0].annotate(str(int(np.round(area_select_ohtsia2_change[month])))+'%',(area_select_ohtsia2_start[month]+0.07,area_select_ohtsia2_end[month]-0.1),color='orange',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_ohtsiv1_change[month])))+'%',(area_select_ohtsiv1_start[month],area_select_ohtsiv1_end[month]+0.1),color='gray',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_ohtsiv2_change[month])))+'%',(area_select_ohtsiv2_start[month]-0.07,area_select_ohtsiv2_end[month]+0.15),color='gray',fontsize=14)
    ax[0,0].annotate(str(int(np.round(area_select_members_change[month])))+'%',(area_select_members_start[month]-0.3,area_select_members_end[month]-0.15),color='purple',fontsize=14)
ax[0,0].set_xlabel('March sea-ice area 2015-2019 (10$^6$ km$^2$)',fontsize=18)
ax[0,0].set_ylabel('March sea-ice area \n 2096-2100 (10$^6$ km$^2$)',fontsize=18)
#ax[0,0].set_xlabel('March sea-ice area 2015-2024 (10$^6$ km$^2$)',fontsize=18)
#ax[0,0].set_ylabel('March sea-ice area \n 2091-2100 (10$^6$ km$^2$)',fontsize=18)
ax[0,0].tick_params(axis='both',labelsize=16)
ax[0,0].grid(linestyle='--')
ax[0,0].set_title('a',loc='left',fontsize=25,fontweight='bold')

# Change in September sea-ice area (2096-2100 vs 2015-2019)
month = 8
ax[0,1].plot(area_mmm_start[month],area_mmm_end[month],'o',color='black',markersize=16)
ax[0,1].errorbar(area_mmm_start[month],area_mmm_end[month],xerr=sd_area_mmm_start[month],yerr=sd_area_mmm_end[month],fmt='ko',capsize=5)
ax[0,1].plot(area_select_sia_15_start[month],area_select_sia_15_end[month],'o',color='blue',markersize=14)
ax[0,1].plot(area_select_siv_15_start[month],area_select_siv_15_end[month],'o',color='green',markersize=14)
ax[0,1].plot(area_select_varsia_15_start[month],area_select_varsia_15_end[month],'X',color='blue',markersize=14)
ax[0,1].plot(area_select_varsiv_15_start[month],area_select_varsiv_15_end[month],'X',color='green',markersize=14)
ax[0,1].plot(area_select_tsia_15_start[month],area_select_tsia_15_end[month],'P',color='blue',markersize=14)
ax[0,1].plot(area_select_tsiv_15_start[month],area_select_tsiv_15_end[month],'P',color='green',markersize=14)
ax[0,1].plot(area_select_aoht_8_start[month],area_select_aoht_8_end[month],'o',color='red',markersize=14)
ax[0,1].plot(area_select_apoht_8_start[month],area_select_apoht_8_end[month],'X',color='red',markersize=14)
ax[0,1].plot(area_select_ohtsia1_start[month],area_select_ohtsia1_end[month],'o',color='orange',markersize=14)
ax[0,1].plot(area_select_ohtsia2_start[month],area_select_ohtsia2_end[month],'X',color='orange',markersize=14)
ax[0,1].plot(area_select_ohtsiv1_start[month],area_select_ohtsiv1_end[month],'o',color='gray',markersize=14)
ax[0,1].plot(area_select_ohtsiv2_start[month],area_select_ohtsiv2_end[month],'X',color='gray',markersize=14)
ax[0,1].plot(area_select_members_start[month],area_select_members_end[month],'o',color='purple',markersize=14)
ax[0,1].plot(area_select_random_start[month],area_select_random_end[month],'o',color='cyan',markersize=10)
ax[0,1].errorbar(area_select_random_start[month],area_select_random_end[month],xerr=sd_area_select_random_start[month],yerr=sd_area_select_random_end[month],fmt='o',color='cyan',capsize=5)
ax[0,1].annotate(str(int(np.round(area_mmm_change[month])))+'%',(area_mmm_start[month]+0.05,area_mmm_end[month]+0.1),color='black',fontsize=16)
if experiment == 'ssp585':
    ax[0,1].annotate(str(int(np.round(area_select_sia_15_change[month])))+'%',(area_select_sia_15_start[month]+0.1,area_select_sia_15_end[month]),color='blue',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_siv_15_change[month])))+'%',(area_select_siv_15_start[month]-0.1,area_select_sia_15_end[month]+0.15),color='green',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_varsia_15_change[month])))+'%',(area_select_varsia_15_start[month]+0.05,area_select_varsia_15_end[month]+0.05),color='blue',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_varsiv_15_change[month])))+'%',(area_select_varsiv_15_start[month]-0.1,area_select_varsiv_15_end[month]+0.08),color='green',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_tsiv_15_change[month])))+'%',(area_select_tsiv_15_start[month]-0.3,area_select_tsiv_15_end[month]+0.1),color='green',fontsize=14)
else:
    ax[0,1].annotate(str(int(np.round(area_select_sia_15_change[month])))+'%',(area_select_sia_15_start[month]-0.2,area_select_sia_15_end[month]-0.3),color='blue',fontsize=13)
    ax[0,1].annotate(str(int(np.round(area_select_siv_15_change[month])))+'%',(area_select_siv_15_start[month]-0.2,area_select_siv_15_end[month]+0.15),color='green',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_varsia_15_change[month])))+'%',(area_select_varsia_15_start[month]+0.05,area_select_varsia_15_end[month]+0.1),color='blue',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_varsiv_15_change[month])))+'%',(area_select_varsiv_15_start[month]+0.05,area_select_varsiv_15_end[month]+0.1),color='green',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_tsia_15_change[month])))+'%',(area_select_tsia_15_start[month]+0.1,area_select_tsia_15_end[month]),color='blue',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_tsiv_15_change[month])))+'%',(area_select_tsiv_15_start[month]-0.1,area_select_tsiv_15_end[month]+0.1),color='green',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_aoht_8_change[month])))+'%',(area_select_aoht_8_start[month],area_select_aoht_8_end[month]-0.3),color='red',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_apoht_8_change[month])))+'%',(area_select_apoht_8_start[month]-0.15,area_select_apoht_8_end[month]+0.15),color='red',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_ohtsia1_change[month])))+'%',(area_select_ohtsia1_start[month]+0.1,area_select_ohtsia1_end[month]),color='orange',fontsize=13)
    ax[0,1].annotate(str(int(np.round(area_select_ohtsia2_change[month])))+'%',(area_select_ohtsia2_start[month]-0.15,area_select_ohtsia2_end[month]-0.3),color='orange',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_ohtsiv1_change[month])))+'%',(area_select_ohtsiv1_start[month]-0.4,area_select_ohtsiv1_end[month]-0.1),color='gray',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_ohtsiv2_change[month])))+'%',(area_select_ohtsiv2_start[month]-0.2,area_select_ohtsiv2_end[month]-0.3),color='gray',fontsize=14)
    ax[0,1].annotate(str(int(np.round(area_select_members_change[month])))+'%',(area_select_members_start[month]-0.3,area_select_members_end[month]-0.2),color='purple',fontsize=12)
ax[0,1].axvline(x=area_obs_start[month],color='black',linestyle='--',linewidth=2)
ax[0,1].set_xlabel('September sea-ice area 2015-2019 (10$^6$ km$^2$)',fontsize=18)
ax[0,1].set_ylabel('September sea-ice area \n 2096-2100 (10$^6$ km$^2$)',fontsize=18)
#ax[0,1].set_xlabel('September sea-ice area 2015-2024 (10$^6$ km$^2$)',fontsize=18)
#ax[0,1].set_ylabel('September sea-ice area \n 2091-2100 (10$^6$ km$^2$)',fontsize=18)
ax[0,1].tick_params(axis='both',labelsize=16)
ax[0,1].grid(linestyle='--')
ax[0,1].set_title('b',loc='left',fontsize=25,fontweight='bold')

# Change in March sea-ice volume (2096-2100 vs 2015-2019)
month = 2
ax[1,0].plot(volume_mmm_start[month],volume_mmm_end[month],'o',color='black',markersize=16)
ax[1,0].errorbar(volume_mmm_start[month],volume_mmm_end[month],xerr=sd_volume_mmm_start[month],yerr=sd_volume_mmm_end[month],fmt='ko',capsize=5)
ax[1,0].plot(volume_select_sia_15_start[month],volume_select_sia_15_end[month],'o',color='blue',markersize=14)
ax[1,0].plot(volume_select_siv_15_start[month],volume_select_siv_15_end[month],'o',color='green',markersize=14)
ax[1,0].plot(volume_select_varsia_15_start[month],volume_select_varsia_15_end[month],'X',color='blue',markersize=14)
ax[1,0].plot(volume_select_varsiv_15_start[month],volume_select_varsiv_15_end[month],'X',color='green',markersize=14)
ax[1,0].plot(volume_select_tsia_15_start[month],volume_select_tsia_15_end[month],'P',color='blue',markersize=14)
ax[1,0].plot(volume_select_tsiv_15_start[month],volume_select_tsiv_15_end[month],'P',color='green',markersize=14)
ax[1,0].plot(volume_select_aoht_8_start[month],volume_select_aoht_8_end[month],'o',color='red',markersize=14)
ax[1,0].plot(volume_select_apoht_8_start[month],volume_select_apoht_8_end[month],'X',color='red',markersize=14)
ax[1,0].plot(volume_select_ohtsia1_start[month],volume_select_ohtsia1_end[month],'o',color='orange',markersize=14)
ax[1,0].plot(volume_select_ohtsia2_start[month],volume_select_ohtsia2_end[month],'X',color='orange',markersize=14)
ax[1,0].plot(volume_select_ohtsiv1_start[month],volume_select_ohtsiv1_end[month],'o',color='gray',markersize=14)
ax[1,0].plot(volume_select_ohtsiv2_start[month],volume_select_ohtsiv2_end[month],'X',color='gray',markersize=14)
ax[1,0].plot(volume_select_members_start[month],volume_select_members_end[month],'o',color='purple',markersize=14)
ax[1,0].plot(volume_select_random_start[month],volume_select_random_end[month],'o',color='cyan',markersize=10)
ax[1,0].errorbar(volume_select_random_start[month],volume_select_random_end[month],xerr=sd_volume_select_random_start[month],yerr=sd_volume_select_random_end[month],fmt='o',color='cyan',capsize=5)
ax[1,0].annotate(str(int(np.round(volume_mmm_change[month])))+'%',(volume_mmm_start[month]+0.1,volume_mmm_end[month]-0.5),color='black',fontsize=16)
if experiment == 'ssp585':
    ax[1,0].annotate(str(int(np.round(volume_select_sia_15_change[month])))+'%',(volume_select_sia_15_start[month]+0.2,volume_select_sia_15_end[month]),color='blue',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_siv_15_change[month])))+'%',(volume_select_siv_15_start[month]-1.,volume_select_siv_15_end[month]+0.1),color='green',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_varsia_15_change[month])))+'%',(volume_select_varsia_15_start[month]+0.05,volume_select_varsia_15_end[month]+0.2),color='blue',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_varsiv_15_change[month])))+'%',(volume_select_varsiv_15_start[month]+0.15,volume_select_varsiv_15_end[month]+0.1),color='green',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_tsia_15_change[month])))+'%',(volume_select_tsia_15_start[month]-0.8,volume_select_tsia_15_end[month]+0.15),color='blue',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_tsiv_15_change[month])))+'%',(volume_select_tsiv_15_start[month]-0.8,volume_select_tsiv_15_end[month]+0.2),color='green',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_aoht_8_change[month])))+'%',(volume_select_aoht_8_start[month]+0.15,volume_select_aoht_8_end[month]-0.3),color='red',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_apoht_8_change[month])))+'%',(volume_select_apoht_8_start[month]-0.15,volume_select_apoht_8_end[month]+0.25),color='red',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_ohtsia1_change[month])))+'%',(volume_select_ohtsia1_start[month]+0.15,volume_select_ohtsia1_end[month]+0.15),color='orange',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_ohtsia2_change[month])))+'%',(volume_select_ohtsia2_start[month]-1.,volume_select_ohtsia2_end[month]),color='orange',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_ohtsiv1_change[month])))+'%',(volume_select_ohtsiv1_start[month]+0.2,volume_select_ohtsiv1_end[month]-0.3),color='gray',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_ohtsiv2_change[month])))+'%',(volume_select_ohtsiv2_start[month]+0.2,volume_select_ohtsiv2_end[month]-0.2),color='gray',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_members_change[month])))+'%',(volume_select_members_start[month]+0.1,volume_select_members_end[month]-0.5),color='purple',fontsize=14)
else:
    ax[1,0].annotate(str(int(np.round(volume_select_sia_15_change[month])))+'%',(volume_select_sia_15_start[month]+0.2,volume_select_sia_15_end[month]-0.3),color='blue',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_siv_15_change[month])))+'%',(volume_select_siv_15_start[month]-1.,volume_select_siv_15_end[month]+0.1),color='green',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_varsia_15_change[month])))+'%',(volume_select_varsia_15_start[month]+0.05,volume_select_varsia_15_end[month]+0.2),color='blue',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_varsiv_15_change[month])))+'%',(volume_select_varsiv_15_start[month]+0.15,volume_select_varsiv_15_end[month]+0.1),color='green',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_tsia_15_change[month])))+'%',(volume_select_tsia_15_start[month]-0.8,volume_select_tsia_15_end[month]+0.15),color='blue',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_tsiv_15_change[month])))+'%',(volume_select_tsiv_15_start[month]-0.8,volume_select_tsiv_15_end[month]+0.2),color='green',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_aoht_8_change[month])))+'%',(volume_select_aoht_8_start[month]+0.15,volume_select_aoht_8_end[month]-0.3),color='red',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_apoht_8_change[month])))+'%',(volume_select_apoht_8_start[month]+0.15,volume_select_apoht_8_end[month]+0.1),color='red',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_ohtsia1_change[month])))+'%',(volume_select_ohtsia1_start[month]-0.1,volume_select_ohtsia1_end[month]-0.5),color='orange',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_ohtsia2_change[month])))+'%',(volume_select_ohtsia2_start[month]-1.,volume_select_ohtsia2_end[month]),color='orange',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_ohtsiv1_change[month])))+'%',(volume_select_ohtsiv1_start[month]-0.5,volume_select_ohtsiv1_end[month]+0.2),color='gray',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_ohtsiv2_change[month])))+'%',(volume_select_ohtsiv2_start[month]-0.3,volume_select_ohtsiv2_end[month]-0.5),color='gray',fontsize=14)
    ax[1,0].annotate(str(int(np.round(volume_select_members_change[month])))+'%',(volume_select_members_start[month]+0.1,volume_select_members_end[month]-0.5),color='purple',fontsize=14)
ax[1,0].axvline(x=volume_piomas_start[month],color='black',linestyle='--',linewidth=2)
ax[1,0].set_xlabel('March sea-ice volume 2015-2019 (10$^3$ km$^3$)',fontsize=18)
ax[1,0].set_ylabel('March sea-ice volume \n 2096-2100 (10$^3$ km$^3$)',fontsize=18)
#ax[1,0].set_xlabel('March sea-ice volume 2015-2024 (10$^3$ km$^3$)',fontsize=18)
#ax[1,0].set_ylabel('March sea-ice volume \n 2091-2100 (10$^3$ km$^3$)',fontsize=18)
ax[1,0].tick_params(axis='both',labelsize=16)
ax[1,0].grid(linestyle='--')
ax[1,0].set_title('c',loc='left',fontsize=25,fontweight='bold')

# Change in September sea-ice volume (2096-2100 vs 2015-2019)
month = 8
if experiment == 'ssp126':
    ax[1,1].plot(volume_mmm_start[month],volume_mmm_end[month],'o',color='black',label='Without selection (32)',markersize=16)
elif experiment == 'ssp585':
    ax[1,1].plot(volume_mmm_start[month],volume_mmm_end[month],'o',color='black',label='Without selection (33)',markersize=16)
ax[1,1].errorbar(volume_mmm_start[month],volume_mmm_end[month],xerr=sd_volume_mmm_start[month],yerr=sd_volume_mmm_end[month],fmt='ko',capsize=5)
ax[1,1].plot(volume_select_sia_15_start[month],volume_select_sia_15_end[month],'o',color='blue',label='Mean sea-ice area (15)',markersize=14)
ax[1,1].plot(volume_select_siv_15_start[month],volume_select_siv_15_end[month],'o',color='green',label='Mean sea-ice volume (15)',markersize=14)
ax[1,1].plot(volume_select_varsia_15_start[month],volume_select_varsia_15_end[month],'X',color='blue',label='Sea-ice area variability (15)',markersize=14)
ax[1,1].plot(volume_select_varsiv_15_start[month],volume_select_varsiv_15_end[month],'X',color='green',label='Sea-ice volume variability (15)',markersize=14)
ax[1,1].plot(volume_select_tsia_15_start[month],volume_select_tsia_15_end[month],'P',color='blue',label='Trend in sea-ice area (15)',markersize=14)
ax[1,1].plot(volume_select_tsiv_15_start[month],volume_select_tsiv_15_end[month],'P',color='green',label='Trend in sea-ice volume (15)',markersize=14)
ax[1,1].plot(volume_select_aoht_8_start[month],volume_select_aoht_8_end[month],'o',color='red',label='Atlantic OHT (8)',markersize=14)
ax[1,1].plot(volume_select_apoht_8_start[month],volume_select_apoht_8_end[month],'X',color='red',label='Atlantic/Pacific OHT (8)',markersize=14)
ax[1,1].plot(volume_select_ohtsia1_start[month],volume_select_ohtsia1_end[month],'o',color='orange',label='Atlantic OHT + sea-ice area (6)',markersize=14)
ax[1,1].plot(volume_select_ohtsia2_start[month],volume_select_ohtsia2_end[month],'X',color='orange',label='Atl/Pac OHT + sea-ice area (5)',markersize=14)
ax[1,1].plot(volume_select_ohtsiv1_start[month],volume_select_ohtsiv1_end[month],'o',color='gray',label='Atlantic OHT + sea-ice volume (3)',markersize=14)
ax[1,1].plot(volume_select_ohtsiv2_start[month],volume_select_ohtsiv2_end[month],'X',color='gray',label='Atl/Pac OHT + sea-ice volume (5)',markersize=14)
ax[1,1].plot(volume_select_members_start[month],volume_select_members_end[month],'o',color='purple',label='>= 5 members (10)',markersize=14)
ax[1,1].plot(volume_select_random_start[month],volume_select_random_end[month],'o',color='cyan',label='Random selection (10)',markersize=10)
ax[1,1].errorbar(volume_select_random_start[month],volume_select_random_end[month],xerr=sd_volume_select_random_start[month],yerr=sd_volume_select_random_end[month],fmt='o',color='cyan',capsize=5)
ax[1,1].annotate(str(int(np.round(volume_mmm_change[month])))+'%',(volume_mmm_start[month]-1.,volume_mmm_end[month]+0.05),color='black',fontsize=16)
if experiment == 'ssp585':
    ax[1,1].annotate(str(int(np.round(volume_select_siv_15_change[month])))+'%',(volume_select_siv_15_start[month]+0.2,volume_select_sia_15_end[month]+0.05),color='green',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_varsia_15_change[month])))+'%',(volume_select_varsia_15_start[month]+0.3,volume_select_varsia_15_end[month]),color='blue',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_varsiv_15_change[month])))+'%',(volume_select_varsiv_15_start[month]+0.3,volume_select_varsiv_15_end[month]),color='green',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_tsiv_15_change[month])))+'%',(volume_select_tsiv_15_start[month],volume_select_tsiv_15_end[month]+0.03),color='green',fontsize=14)
else:
    ax[1,1].annotate(str(int(np.round(volume_select_sia_15_change[month])))+'%',(volume_select_sia_15_start[month]+0.2,volume_select_sia_15_end[month]),color='blue',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_siv_15_change[month])))+'%',(volume_select_siv_15_start[month]-0.5,volume_select_siv_15_end[month]+0.3),color='green',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_varsia_15_change[month])))+'%',(volume_select_varsia_15_start[month]+0.05,volume_select_varsia_15_end[month]+0.2),color='blue',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_varsiv_15_change[month])))+'%',(volume_select_varsiv_15_start[month]+0.15,volume_select_varsiv_15_end[month]+0.1),color='green',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_tsia_15_change[month])))+'%',(volume_select_tsia_15_start[month]-0.8,volume_select_tsia_15_end[month]+0.15),color='blue',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_tsiv_15_change[month])))+'%',(volume_select_tsiv_15_start[month]-0.8,volume_select_tsiv_15_end[month]+0.2),color='green',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_aoht_8_change[month])))+'%',(volume_select_aoht_8_start[month],volume_select_aoht_8_end[month]-0.3),color='red',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_apoht_8_change[month])))+'%',(volume_select_apoht_8_start[month]+0.15,volume_select_apoht_8_end[month]+0.1),color='red',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_ohtsia1_change[month])))+'%',(volume_select_ohtsia1_start[month]-0.4,volume_select_ohtsia1_end[month]+0.15),color='orange',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_ohtsia2_change[month])))+'%',(volume_select_ohtsia2_start[month],volume_select_ohtsia2_end[month]-0.4),color='orange',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_ohtsiv1_change[month])))+'%',(volume_select_ohtsiv1_start[month]+0.2,volume_select_ohtsiv1_end[month]-0.3),color='gray',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_ohtsiv2_change[month])))+'%',(volume_select_ohtsiv2_start[month]-0.5,volume_select_ohtsiv2_end[month]-0.4),color='gray',fontsize=14)
    ax[1,1].annotate(str(int(np.round(volume_select_members_change[month])))+'%',(volume_select_members_start[month]+0.15,volume_select_members_end[month]+0.15),color='purple',fontsize=14)
ax[1,1].axvline(x=volume_piomas_start[month],color='black',label='Obs./Reanalysis',linestyle='--',linewidth=2)
ax[1,1].legend(shadow=True,frameon=False,fontsize=16,bbox_to_anchor=(1.05,-0.18),ncol=3)
ax[1,1].set_xlabel('September sea-ice volume 2015-2019 (10$^3$ km$^3$)',fontsize=18)
ax[1,1].set_ylabel('September sea-ice volume \n 2096-2100 (10$^3$ km$^3$)',fontsize=18)
#ax[1,1].set_xlabel('September sea-ice volume 2015-2024 (10$^3$ km$^3$)',fontsize=18)
#ax[1,1].set_ylabel('September sea-ice volume \n 2091-2100 (10$^3$ km$^3$)',fontsize=18)
ax[1,1].tick_params(axis='both',labelsize=16)
ax[1,1].grid(linestyle='--')
ax[1,1].set_title('d',loc='left',fontsize=25,fontweight='bold')

# Save figure
if save_fig == True:
    if experiment == 'ssp585':
        fig.savefig(dir_output + 'fig3.png')
    else:
        fig.savefig(dir_output + 'figS5.png')
