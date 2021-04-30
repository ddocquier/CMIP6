#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 4: Plot date of first ice-free Arctic in September for each model selection criterion based on SSP5-8.5
    Fig. S6: Same as Fig. 4 for SSP1-2.6
PROGRAMMER
    D. Docquier
LAST UPDATEs
    30/04/2021
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt

# Option
experiment = 'ssp585' # ssp126; ssp585
save_fig = True

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/' + str(experiment) + '/SIA-SIV/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6_Paper/'

# Time parameters
nmy = 12

# Function to compute date of first ice-free Arctic
def compute_date(years,var):
    var_sep = var[:,8]
    year = np.nan
    y = 0
    while y < nyears:
        if var_sep[y] < 1.:
            year = y + 2015
            break
        y = y + 1
    return year

# Load multi-model means from model selections
if experiment == 'ssp585':
    filename = dir_output + 'area_mmm_ssp585.npy'
else:
    filename = dir_output + 'area_mmm_ssp126.npy'
area_mmm,area_select_sia,area_select_siv,area_select_varsia,area_select_varsiv, \
    area_select_tsia,area_select_tsiv,area_select_aoht,area_select_apoht,area_select_ohtsia1, \
    area_select_ohtsia2,area_select_ohtsiv1,area_select_ohtsiv2,area_select_members = np.load(filename,allow_pickle=True)
nm = np.size(area_mmm)
nyears = int(nm / nmy)
mmm_date = compute_date(nyears,area_mmm)
sia_date = compute_date(nyears,area_select_sia)
siv_date = compute_date(nyears,area_select_siv)
varsia_date = compute_date(nyears,area_select_varsia)
varsiv_date = compute_date(nyears,area_select_varsiv)
tsia_date = compute_date(nyears,area_select_tsia)
tsiv_date = compute_date(nyears,area_select_tsiv)
aoht_date = compute_date(nyears,area_select_aoht)
apoht_date = compute_date(nyears,area_select_apoht)
ohtsia1_date = compute_date(nyears,area_select_ohtsia1)
ohtsia2_date = compute_date(nyears,area_select_ohtsia2)
ohtsiv1_date = compute_date(nyears,area_select_ohtsiv1)
ohtsiv2_date = compute_date(nyears,area_select_ohtsiv2)
members_date = compute_date(nyears,area_select_members)

# Load sea-ice area and volume AWI-CM-1-1-MR
filename = dir_input + 'SIA-SIV_SImon_AWI-CM-1-1-MR_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_awicm,volume_awicm = np.load(filename)
area_awicm_date = compute_date(nyears,area_awicm)

# Load sea-ice area and volume BCC-CSM2-MR
filename = dir_input + 'SIA-SIV_SImon_BCC-CSM2-MR_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_bcccsm2mr,volume_bcccsm2mr = np.load(filename)
area_bcccsm2mr_date = compute_date(nyears,area_bcccsm2mr)

# Load sea-ice area and volume CAMS-CSM1-0
filename = dir_input + 'SIA-SIV_SImon_CAMS-CSM1-0_' + str(experiment) + '_ensmean_gn_2015-2099.npy'
area_camscsm_init,volume_camscsm_init,notused,notused = np.load(filename)
area_camscsm = np.zeros((nyears,12)) * np.nan
volume_camscsm = np.zeros((nyears,12)) * np.nan
area_camscsm[0:85,:] = area_camscsm_init
volume_camscsm[0:85,:] = volume_camscsm_init
area_camscsm_date = compute_date(nyears,area_camscsm)

# Load sea-ice area and volume FGOALS-f3-L
filename = dir_input + 'SIA-SIV_SImon_FGOALS-f3-L_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_fgoalsf3l,volume_fgoalsf3l = np.load(filename)
area_fgoalsf3l_date = compute_date(nyears,area_fgoalsf3l)

# Load sea-ice area and volume FGOALS-g3
filename = dir_input + 'SIA-SIV_SImon_FGOALS-g3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_fgoalsg3,volume_fgoalsg3,notused,notused = np.load(filename)
area_fgoalsg3_date = compute_date(nyears,area_fgoalsg3)

# Load sea-ice area and volume CanESM5
filename = dir_input + 'SIA-SIV_SImon_CanESM5_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_canesm5,volume_canesm5,notused,notused = np.load(filename)
area_canesm5_date = compute_date(nyears,area_canesm5)

# Load sea-ice area and volume CanESM5-CanOE
filename = dir_input + 'SIA-SIV_SImon_CanESM5-CanOE_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_canesm5canoe,volume_canesm5canoe,notused,notused = np.load(filename)
area_canesm5canoe_date = compute_date(nyears,area_canesm5canoe)

# Load sea-ice area and volume CNRM-CM6-1
filename = dir_input + 'SIA-SIV_SImon_CNRM-CM6-1_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cnrmcm6,volume_cnrmcm6,notused,notused = np.load(filename)
area_cnrmcm6_date = compute_date(nyears,area_cnrmcm6)

# Load sea-ice area and volume CNRM-CM6-1-HR
filename = dir_input + 'SIA-SIV_SImon_CNRM-CM6-1-HR_' + str(experiment) + '_r1i1p1f2_gn_2015-2100.npy'
area_cnrmcm6hr,volume_cnrmcm6hr = np.load(filename)
area_cnrmcm6hr_date = compute_date(nyears,area_cnrmcm6hr)

# Load sea-ice area and volume ACCESS-ESM1-5
filename = dir_input + 'SIA-SIV_SImon_ACCESS-ESM1-5_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_accessesm,volume_accessesm,notused,notused = np.load(filename)
area_accessesm_date = compute_date(nyears,area_accessesm)

# Load sea-ice area and volume ACCESS-CM2
filename = dir_input + 'SIA-SIV_SImon_ACCESS-CM2_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_accesscm,volume_accesscm,notused,notused = np.load(filename)
area_accesscm_date = compute_date(nyears,area_accesscm)

# Load sea-ice area and volume MPI-ESM1-2-LR
filename = dir_input + 'SIA-SIV_SImon_MPI-ESM1-2-LR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_mpiesmlr,volume_mpiesmlr,notused,notused = np.load(filename)
area_mpiesmlr_date = compute_date(nyears,area_mpiesmlr)

# Load sea-ice area and volume MPI-ESM1-2-HR
filename = dir_input + 'SIA-SIV_SImon_MPI-ESM1-2-HR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_mpiesmhr,volume_mpiesmhr,notused,notused = np.load(filename)
area_mpiesmhr_date = compute_date(nyears,area_mpiesmhr)

# Load sea-ice area and volume EC-Earth3
filename = dir_input + 'SIA-SIV_SImon_EC-Earth3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ecearth3,volume_ecearth3,notused,notused = np.load(filename)
area_ecearth3_date = compute_date(nyears,area_ecearth3)

# Load sea-ice area and volume EC-Earth3-Veg
filename = dir_input + 'SIA-SIV_SImon_EC-Earth3-Veg_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ecearth3veg,volume_ecearth3veg,notused,notused = np.load(filename)
area_ecearth3veg_date = compute_date(nyears,area_ecearth3veg)

# Load sea-ice area and volume FIO-ESM-2-0
filename = dir_input + 'SIA-SIV_SImon_FIO-ESM-2-0_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_fioesm,volume_fioesm,notused,notused = np.load(filename)
area_fioesm_date = compute_date(nyears,area_fioesm)

# Load sea-ice area and volume INM-CM4-8
filename = dir_input + 'SIA-SIV_SImon_INM-CM4-8_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_inmcm48,volume_inmcm48 = np.load(filename)
area_inmcm48_date = compute_date(nyears,area_inmcm48)

# Load sea-ice area and volume INM-CM5-0
filename = dir_input + 'SIA-SIV_SImon_INM-CM5-0_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_inmcm50,volume_inmcm50 = np.load(filename)
area_inmcm50_date = compute_date(nyears,area_inmcm50)

# Load sea-ice area and volume IPSL-CM6A-LR
filename = dir_input + 'SIA-SIV_SImon_IPSL-CM6A-LR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ipslcm6alr,volume_ipslcm6alr,notused,notused = np.load(filename)
area_ipslcm6alr_date = compute_date(nyears,area_ipslcm6alr)

# Load sea-ice area and volume KIOST-ESM
filename = dir_input + 'SIA-SIV_SImon_KIOST-ESM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_kiostesm,volume_kiostesm = np.load(filename)
area_kiostesm_date = compute_date(nyears,area_kiostesm)

# Load sea-ice area and volume MIROC6
filename = dir_input + 'SIA-SIV_SImon_MIROC6_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_miroc6,volume_miroc6,notused,notused = np.load(filename)
area_miroc6_date = compute_date(nyears,area_miroc6)

# Load sea-ice area and volume MIROC-ES2L
filename = dir_input + 'SIA-SIV_SImon_MIROC-ES2L_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_miroces2l,volume_miroces2l,notused,notused = np.load(filename)
area_miroces2l_date = compute_date(nyears,area_miroces2l)

# Load sea-ice area and volume HadGEM3-GC31-LL
filename = dir_input + 'SIA-SIV_SImon_HadGEM3-GC31-LL_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_hadgem3ll,volume_hadgem3ll,notused,notused = np.load(filename)
area_hadgem3ll_date = compute_date(nyears,area_hadgem3ll)

# Load sea-ice area and volume HadGEM3-GC31-MM
filename = dir_input + 'SIA-SIV_SImon_HadGEM3-GC31-MM_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_hadgem3mm,volume_hadgem3mm,notused,notused = np.load(filename)
area_hadgem3mm_date = compute_date(nyears,area_hadgem3mm)

# Load sea-ice area and volume UKESM1-0-LL
filename = dir_input + 'SIA-SIV_SImon_UKESM1-0-LL_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ukesmll,volume_ukesmll,notused,notused = np.load(filename)
area_ukesmll_date = compute_date(nyears,area_ukesmll)

# Load sea-ice area and volume MRI-ESM2-0
filename = dir_input + 'SIA-SIV_SImon_MRI-ESM2-0_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_mriesm,volume_mriesm = np.load(filename)
area_mriesm_date = compute_date(nyears,area_mriesm)

# Load sea-ice area and volume CESM2
filename = dir_input + 'SIA-SIV_SImon_CESM2_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cesm2,volume_cesm2,notused,notused = np.load(filename)
area_cesm2_date = compute_date(nyears,area_cesm2)

# Load sea-ice area and volume CESM2-WACCM
filename = dir_input + 'SIA-SIV_SImon_CESM2-WACCM_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cesm2waccm,volume_cesm2waccm,notused,notused = np.load(filename)
area_cesm2waccm_date = compute_date(nyears,area_cesm2waccm)

# Load sea-ice area and volume NorESM2-LM
filename = dir_input + 'SIA-SIV_SImon_NorESM2-LM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_noresm2lm,volume_noresm2lm = np.load(filename)
area_noresm2lm_date = compute_date(nyears,area_noresm2lm)

# Load sea-ice area and volume NorESM2-MM
filename = dir_input + 'SIA-SIV_SImon_NorESM2-MM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_noresm2mm,volume_noresm2mm = np.load(filename)
area_noresm2mm_date = compute_date(nyears,area_noresm2mm)

# Load sea-ice area and volume GFDL-CM4
if experiment == 'ssp585':
    filename = dir_input + 'SIA-SIV_SImon_GFDL-CM4_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
    area_gfdlcm4,volume_gfdlcm4 = np.load(filename)
    area_gfdlcm4_date = compute_date(nyears,area_gfdlcm4)

# Load sea-ice area and volume GFDL-ESM4
filename = dir_input + 'SIA-SIV_SImon_GFDL-ESM4_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_gfdlesm4,volume_gfdlesm4 = np.load(filename)
area_gfdlesm4_date = compute_date(nyears,area_gfdlesm4)

# Load sea-ice area and volume NESM3
filename = dir_input + 'SIA-SIV_SImon_NESM3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_nesm3,volume_nesm3,notused,notused = np.load(filename)
area_nesm3_date = compute_date(nyears,area_nesm3)

# All models
if experiment == 'ssp585':
    array_mmm = [area_awicm_date,area_bcccsm2mr_date,area_camscsm_date,area_fgoalsf3l_date,area_fgoalsg3_date,area_canesm5_date,area_canesm5canoe_date,area_cnrmcm6_date,area_cnrmcm6hr_date,area_accessesm_date,area_accesscm_date,area_mpiesmlr_date,area_mpiesmhr_date,area_ecearth3_date,area_ecearth3veg_date,area_fioesm_date,area_inmcm48_date,area_inmcm50_date,area_ipslcm6alr_date,area_kiostesm_date,area_miroc6_date,area_miroces2l_date,area_hadgem3ll_date,area_hadgem3mm_date,area_ukesmll_date,area_mriesm_date,area_cesm2_date,area_cesm2waccm_date,area_noresm2lm_date,area_noresm2mm_date,area_gfdlcm4_date,area_gfdlesm4_date,area_nesm3_date]
else:
    array_mmm = [area_awicm_date,area_bcccsm2mr_date,area_camscsm_date,area_fgoalsf3l_date,area_fgoalsg3_date,area_canesm5_date,area_canesm5canoe_date,area_cnrmcm6_date,area_cnrmcm6hr_date,area_accessesm_date,area_accesscm_date,area_mpiesmlr_date,area_mpiesmhr_date,area_ecearth3_date,area_ecearth3veg_date,area_fioesm_date,area_inmcm48_date,area_inmcm50_date,area_ipslcm6alr_date,area_kiostesm_date,area_miroc6_date,area_miroces2l_date,area_hadgem3ll_date,area_hadgem3mm_date,area_ukesmll_date,area_mriesm_date,area_cesm2_date,area_cesm2waccm_date,area_noresm2lm_date,area_noresm2mm_date,area_gfdlesm4_date,area_nesm3_date]
#mmm_date = np.nanmean(array_mmm)
sd_mmm_date = np.nanstd(array_mmm)
print('MMM:',array_mmm)

# Good mean SIA (15 best models)
if experiment == 'ssp585':
    array_sia = [area_cesm2waccm_date,area_noresm2lm_date,area_accessesm_date,area_gfdlesm4_date,area_mpiesmlr_date,area_cnrmcm6hr_date,area_hadgem3mm_date,area_ipslcm6alr_date,area_gfdlcm4_date,area_inmcm50_date,area_cnrmcm6_date,area_ecearth3_date,area_ecearth3veg_date,area_hadgem3ll_date,area_miroces2l_date]
else:
    array_sia = [area_cesm2waccm_date,area_noresm2lm_date,area_accessesm_date,area_gfdlesm4_date,area_mpiesmlr_date,area_cnrmcm6hr_date,area_hadgem3mm_date,area_ipslcm6alr_date,area_inmcm50_date,area_cnrmcm6_date,area_ecearth3_date,area_ecearth3veg_date,area_hadgem3ll_date,area_miroces2l_date,area_mriesm_date]
#sia_date = np.nanmean(array_sia)
sd_sia_date = np.nanstd(array_sia)
print('SIA:',array_sia)

# Good mean SIV (15 best models)
if experiment == 'ssp585':
    array_siv = [area_cesm2waccm_date,area_ipslcm6alr_date,area_fgoalsf3l_date,area_camscsm_date,area_miroc6_date,area_mpiesmhr_date,area_fioesm_date,area_accessesm_date,area_mpiesmlr_date,area_awicm_date,area_cesm2_date,area_gfdlcm4_date,area_hadgem3mm_date,area_mriesm_date,area_accesscm_date]
else:
    array_siv = [area_cesm2waccm_date,area_ipslcm6alr_date,area_fgoalsf3l_date,area_camscsm_date,area_miroc6_date,area_mpiesmhr_date,area_fioesm_date,area_accessesm_date,area_mpiesmlr_date,area_awicm_date,area_cesm2_date,area_hadgem3mm_date,area_mriesm_date,area_accesscm_date,area_noresm2lm_date]
#siv_date = np.nanmean(array_siv)
sd_siv_date = np.nanstd(array_siv)
print('SIV:',array_siv)

# Good trend in SIA (15 best models)
if experiment == 'ssp585':
    array_tsia = [area_hadgem3mm_date,area_cesm2waccm_date,area_fgoalsf3l_date,area_mpiesmhr_date,area_cesm2_date,area_mriesm_date,area_accessesm_date,area_bcccsm2mr_date,area_fioesm_date,area_gfdlesm4_date,area_ipslcm6alr_date,area_mpiesmlr_date,area_cnrmcm6hr_date,area_miroc6_date,area_gfdlcm4_date]
else:
    array_tsia = [area_hadgem3mm_date,area_cesm2waccm_date,area_fgoalsf3l_date,area_mpiesmhr_date,area_cesm2_date,area_mriesm_date,area_accessesm_date,area_bcccsm2mr_date,area_fioesm_date,area_gfdlesm4_date,area_ipslcm6alr_date,area_mpiesmlr_date,area_cnrmcm6hr_date,area_miroc6_date,area_hadgem3ll_date]
#tsia_date = np.nanmean(array_tsia)
sd_tsia_date = np.nanstd(array_tsia)
print('TSIA:',array_tsia)

# Compute multi-model mean - good trend in SIV (15 best models)
if experiment == 'ssp585':
    array_tsiv = [area_ipslcm6alr_date,area_cesm2_date,area_fioesm_date,area_nesm3_date,area_noresm2mm_date,area_noresm2lm_date,area_mpiesmhr_date,area_miroc6_date,area_accesscm_date,area_gfdlcm4_date,area_mpiesmlr_date,area_fgoalsf3l_date,area_gfdlesm4_date,area_kiostesm_date,area_ecearth3_date]
else:
    array_tsiv = [area_ipslcm6alr_date,area_cesm2_date,area_fioesm_date,area_nesm3_date,area_noresm2mm_date,area_noresm2lm_date,area_mpiesmhr_date,area_miroc6_date,area_accesscm_date,area_mpiesmlr_date,area_fgoalsf3l_date,area_gfdlesm4_date,area_kiostesm_date,area_ecearth3_date,area_mriesm_date]
#tsiv_date = np.nanmean(array_tsiv)
sd_tsiv_date = np.nanstd(array_tsiv)
print('TSIV:',array_tsiv)

# Good SIA variability (15 best models)
if experiment == 'ssp585':
    array_varsia = [area_gfdlcm4_date,area_accessesm_date,area_inmcm50_date,area_noresm2mm_date,area_gfdlesm4_date,area_canesm5canoe_date,area_cnrmcm6hr_date,area_accesscm_date,area_hadgem3ll_date,area_nesm3_date,area_fgoalsf3l_date,area_fgoalsg3_date,area_inmcm48_date,area_bcccsm2mr_date,area_ecearth3_date]
else:
    array_varsia = [area_accessesm_date,area_inmcm50_date,area_noresm2mm_date,area_gfdlesm4_date,area_canesm5canoe_date,area_cnrmcm6hr_date,area_accesscm_date,area_hadgem3ll_date,area_nesm3_date,area_fgoalsf3l_date,area_fgoalsg3_date,area_inmcm48_date,area_bcccsm2mr_date,area_ecearth3_date,area_mriesm_date]
#varsia_date = np.nanmean(array_varsia)
sd_varsia_date = np.nanstd(array_varsia)
print('VSIA:',array_varsia)

# Good SIV variability (15 best models)
if experiment == 'ssp585':
    array_varsiv = [area_gfdlcm4_date,area_canesm5_date,area_noresm2lm_date,area_noresm2mm_date,area_fioesm_date,area_hadgem3mm_date,area_cesm2waccm_date,area_miroces2l_date,area_accesscm_date,area_gfdlesm4_date,area_cesm2_date,area_ecearth3veg_date,area_ecearth3_date,area_fgoalsf3l_date,area_camscsm_date]
else:
    array_varsiv = [area_canesm5_date,area_noresm2lm_date,area_noresm2mm_date,area_fioesm_date,area_hadgem3mm_date,area_cesm2waccm_date,area_miroces2l_date,area_accesscm_date,area_gfdlesm4_date,area_cesm2_date,area_ecearth3veg_date,area_ecearth3_date,area_fgoalsf3l_date,area_camscsm_date,area_ukesmll_date]
#varsiv_date = np.nanmean(array_varsiv)
sd_varsiv_date = np.nanstd(array_varsiv)
print('VSIV:',array_varsiv)

# Good AOHT 26N and 57N (8 best models)
array_aoht = [area_hadgem3mm_date,area_ecearth3veg_date,area_mriesm_date,area_mpiesmlr_date,area_ecearth3_date,area_hadgem3ll_date,area_ukesmll_date,area_cnrmcm6_date]
#aoht_date = np.nanmean(array_aoht)
sd_aoht_date = np.nanstd(array_aoht)
print('AOHT:',array_aoht)

# Good AOHT 70N and POHT 60N (8 best models)
array_apoht = [area_mpiesmlr_date,area_ukesmll_date,area_mpiesmhr_date,area_hadgem3ll_date,area_hadgem3mm_date,area_ipslcm6alr_date,area_mriesm_date,area_canesm5_date]
#apoht_date = np.nanmean(array_apoht)
sd_apoht_date = np.nanstd(array_apoht)
print('APOHT:',array_apoht)

# Good AOHT 26N+57N and SIA (among the best 15 SIA models and best 8 OHT models)
array_ohtsia1 = [area_hadgem3mm_date,area_ecearth3veg_date,area_mpiesmlr_date,area_ecearth3_date,area_hadgem3ll_date,area_cnrmcm6_date]
#ohtsia1_date = np.nanmean(array_ohtsia1)
sd_ohtsia1_date = np.nanstd(array_ohtsia1)
print('AOHT/SIA:',array_ohtsia1)

# Good AOHT 70N + POHT 60N and SIA (among the best 15 SIA models and best 8 OHT models)
array_ohtsia2 = [area_mpiesmlr_date,area_hadgem3ll_date,area_hadgem3mm_date,area_ipslcm6alr_date]
#ohtsia2_date = np.nanmean(array_ohtsia2)
sd_ohtsia2_date = np.nanstd(array_ohtsia2)
print('APOHT/SIA:',array_ohtsia2)

# Good AOHT 26N+57N and SIV (among the best 15 SIV models and best 8 OHT models)
array_ohtsiv1 = [area_hadgem3mm_date,area_mriesm_date,area_mpiesmlr_date]
#ohtsiv1_date = np.nanmean(array_ohtsiv1)
sd_ohtsiv1_date = np.nanstd(array_ohtsiv1)
print('AOHT/SIV:',array_ohtsiv1)

# Good AOHT 70N + POHT 60N and SIV (among the best 15 SIV models and best 8 OHT models)
array_ohtsiv2 = [area_mpiesmlr_date,area_mpiesmhr_date,area_hadgem3mm_date,area_ipslcm6alr_date,area_mriesm_date]
#ohtsiv2_date = np.nanmean(array_ohtsiv2)
sd_ohtsiv2_date = np.nanstd(array_ohtsiv2)
print('APOHT/SIV:',array_ohtsiv2)

# Compute multi-model mean - >=5 members
if experiment == 'ssp585':
    array_members = [area_canesm5_date,area_cnrmcm6_date,area_mpiesmlr_date,area_ecearth3_date,area_ecearth3veg_date,area_ipslcm6alr_date,area_miroc6_date,area_ukesmll_date,area_cesm2_date,area_cesm2waccm_date]
else:
    array_members = [area_canesm5_date,area_cnrmcm6_date,area_mpiesmlr_date,area_ecearth3_date,area_ecearth3veg_date,area_ipslcm6alr_date,area_miroc6_date,area_ukesmll_date,area_cesm2_date]
#members_date = np.nanmean(array_members)
sd_members_date = np.nanstd(array_members)
print('Members:',array_members)

# Print
print('WS:',mmm_date)
print('SIA:',sia_date)
print('SIV:',siv_date)
print('VSIA:',varsia_date)
print('VSIV:',varsiv_date)
print('TSIA:',tsia_date)
print('TSIV:',tsiv_date)
print('AOHT:',aoht_date)
print('APOHT:',apoht_date)
print('AOHT+SIA:',ohtsia1_date)
print('APOHT+SIA:',ohtsia2_date)
print('AOHT+SIV:',ohtsiv1_date)
print('APOHT+SIV:',ohtsiv2_date)
print('Members:',members_date)


# Labels
if experiment == 'ssp585':
    name_xticks = ['WS \n 4','C1 \n 0','C2 \n 1','C3 \n 2','C4 \n 1','C5 \n 0','C6 \n 1','C7 \n 0',
               'C8 \n 0','C9 \n 0','C10 \n 0','C11 \n 0','C12 \n 0','C13 \n 0']
else:
    name_xticks = ['WS \n 12','C1 \n 6','C2 \n 4','C3 \n 6','C4 \n 4','C5 \n 4','C6 \n 6','C7 \n 2',
               'C8 \n 1','C9 \n 2','C10 \n 1','C11 \n 1','C12 \n 1','C13 \n 3'] 

# Plot date of first ice-free Arctic in September
index = np.arange(14)
fig,ax = plt.subplots(figsize=(16,12))
fig.subplots_adjust(left=0.1,bottom=0.25,right=0.95,top=0.95,wspace=None,hspace=None)
if experiment == 'ssp126':
    ax.plot(index[0],mmm_date,'o',color='black',label='WS: Without selection (32)',markersize=16)
elif experiment == 'ssp585':
    ax.plot(index[0],mmm_date,'o',color='black',label='WS: Without selection (33)',markersize=16)
ax.plot(index[1],sia_date,'o',color='blue',label='C1: Mean sea-ice area (15)',markersize=14)
ax.plot(index[2],siv_date,'o',color='green',label='C2: Mean sea-ice volume (15)',markersize=14)
ax.plot(index[3],varsia_date,'X',color='blue',label='C3: Sea-ice area variability (15)',markersize=14)
ax.plot(index[4],varsiv_date,'X',color='green',label='C4: Sea-ice volume variability (15)',markersize=14)
ax.plot(index[5],tsia_date,'P',color='blue',label='C5: Trend in sea-ice area (15)',markersize=14)
ax.plot(index[6],tsiv_date,'P',color='green',label='C6: Trend in sea-ice volume (15)',markersize=14)
ax.plot(index[7],aoht_date,'o',color='red',label='C7: Atlantic OHT (8)',markersize=14)
ax.plot(index[8],apoht_date,'X',color='red',label='C8: Atlantic/Pacific OHT (8)',markersize=14)
ax.plot(index[9],ohtsia1_date,'o',color='orange',label='C9: Atlantic OHT + sea-ice area (6)',markersize=14)
ax.plot(index[10],ohtsia2_date,'X',color='orange',label='C10: Atl/Pac OHT + sea-ice area (4)',markersize=14)
ax.plot(index[11],ohtsiv1_date,'o',color='gray',label='C11: Atlantic OHT + sea-ice volume (3)',markersize=14)
ax.plot(index[12],ohtsiv2_date,'X',color='gray',label='C12: Atl/Pac OHT + sea-ice volume (5)',markersize=14)
if experiment == 'ssp585':
    ax.plot(index[13],members_date,'o',color='purple',label='C13: >= 5 members (10)',markersize=14)
elif experiment == 'ssp126':
    ax.plot(index[13],members_date,'o',color='purple',label='C13: >= 5 members (9)',markersize=14)
for i in np.arange(np.size(array_mmm)):
    ax.plot(index[0],array_mmm[i],'.',color='black',markersize=8)
for i in np.arange(np.size(array_sia)):
    ax.plot(index[1],array_sia[i],'.',color='blue',markersize=8)
for i in np.arange(np.size(array_siv)):
    ax.plot(index[2],array_siv[i],'.',color='green',markersize=8)
for i in np.arange(np.size(array_varsia)):
    ax.plot(index[3],array_varsia[i],'.',color='blue',markersize=8)
for i in np.arange(np.size(array_varsiv)):
    ax.plot(index[4],array_varsiv[i],'.',color='green',markersize=8)
for i in np.arange(np.size(array_tsia)):
    ax.plot(index[5],array_tsia[i],'.',color='blue',markersize=8)
for i in np.arange(np.size(array_tsiv)):
    ax.plot(index[6],array_tsiv[i],'.',color='green',markersize=8)
for i in np.arange(np.size(array_aoht)):
    ax.plot(index[7],array_aoht[i],'.',color='red',markersize=8)
for i in np.arange(np.size(array_apoht)):
    ax.plot(index[8],array_apoht[i],'.',color='red',markersize=8)
for i in np.arange(np.size(array_ohtsia1)):
    ax.plot(index[9],array_ohtsia1[i],'.',color='orange',markersize=8)
for i in np.arange(np.size(array_ohtsia2)):
    ax.plot(index[10],array_ohtsia2[i],'.',color='orange',markersize=8)
for i in np.arange(np.size(array_ohtsiv1)):
    ax.plot(index[11],array_ohtsiv1[i],'.',color='gray',markersize=8)
for i in np.arange(np.size(array_ohtsiv2)):
    ax.plot(index[12],array_ohtsiv2[i],'.',color='gray',markersize=8)
for i in np.arange(np.size(array_members)):
    ax.plot(index[13],array_members[i],'.',color='purple',markersize=8)
#ax.errorbar(index[0],mmm_date,yerr=sd_mmm_date,fmt='o',color='black',capsize=5)
#ax.errorbar(index[1],sia_date,yerr=sd_sia_date,fmt='o',color='blue',capsize=5)
#ax.errorbar(index[2],siv_date,yerr=sd_siv_date,fmt='o',color='green',capsize=5)
#ax.errorbar(index[3],varsia_date,yerr=sd_varsia_date,fmt='o',color='blue',capsize=5)
#ax.errorbar(index[4],varsiv_date,yerr=sd_varsiv_date,fmt='o',color='green',capsize=5)
#ax.errorbar(index[5],tsia_date,yerr=sd_tsia_date,fmt='o',color='blue',capsize=5)
#ax.errorbar(index[6],tsiv_date,yerr=sd_tsiv_date,fmt='o',color='green',capsize=5)
#ax.errorbar(index[7],aoht_date,yerr=sd_aoht_date,fmt='o',color='red',capsize=5)
#ax.errorbar(index[8],apoht_date,yerr=sd_apoht_date,fmt='o',color='red',capsize=5)
#ax.errorbar(index[9],ohtsia1_date,yerr=sd_ohtsia1_date,fmt='o',color='orange',capsize=5)
#ax.errorbar(index[10],ohtsia2_date,yerr=sd_ohtsia2_date,fmt='o',color='orange',capsize=5)
#ax.errorbar(index[11],ohtsiv1_date,yerr=sd_ohtsiv1_date,fmt='o',color='gray',capsize=5)
#ax.errorbar(index[12],ohtsiv2_date,yerr=sd_ohtsiv2_date,fmt='o',color='gray',capsize=5)
#ax.errorbar(index[13],members_date,yerr=sd_members_date,fmt='o',color='purple',capsize=5)
ax.legend(shadow=True,frameon=False,fontsize=15,bbox_to_anchor=(1,-0.15),ncol=3)
ax.set_ylabel('Date of first September ice-free Arctic',fontsize=24)
ax.set_xlabel('Model selection',fontsize=24)
ax.set_xticks(index)
ax.set_xticklabels(name_xticks)
ax.set_yticks(np.arange(2020, 2100.1, 10))
ax.tick_params(axis='both',labelsize=20)
ax.set_ylim([2014,2100])
ax.grid(linestyle='--')

# Save figure
if save_fig == True:
    if experiment == 'ssp585':
        fig.savefig(dir_output + 'fig4.png')
    else:
        fig.savefig(dir_output + 'figS6.png')