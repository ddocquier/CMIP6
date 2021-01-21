#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
GOAL
    Fig. 2: Plot Arctic sea-ice area (10^6 km^2) for selected CMIP6 models based on SSP5-8.5
    Fig. S2: Same as Fig. 2 for SSP1-2.6
    Fig. S3: Plot Arctic sea-ice volume (10^3 km^3) for selected CMIP6 models based on SSP5-8.5
    Fig. S4: Same as Fig. S3 for SSP1-2.6
PROGRAMMER
    D. Docquier
LAST UPDATEs
    18/12/2020
'''

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt

# Option
experiment = 'ssp585' # ssp126; ssp585
save_fig = False

# Working directories
dir_input = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6/' + str(experiment) + '/SIA-SIV/'
dir_output = '/nobackup/rossby24/proj/rossby/joint_exp/oseaice/CMIP6_Paper/'

# Time parameters
nmy = 12

# Load sea-ice area and volume AWI-CM-1-1-MR
filename = dir_input + 'SIA-SIV_SImon_AWI-CM-1-1-MR_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_awicm,volume_awicm = np.load(filename)

# Load sea-ice area and volume BCC-CSM2-MR
filename = dir_input + 'SIA-SIV_SImon_BCC-CSM2-MR_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_bcccsm2mr,volume_bcccsm2mr = np.load(filename)
nm = np.size(area_bcccsm2mr)
nyears = int(nm / nmy)

# Load sea-ice area and volume CAMS-CSM1-0
filename = dir_input + 'SIA-SIV_SImon_CAMS-CSM1-0_' + str(experiment) + '_ensmean_gn_2015-2099.npy'
area_camscsm_init,volume_camscsm_init,notused,notused = np.load(filename)
area_camscsm = np.zeros((nyears,12)) * np.nan
volume_camscsm = np.zeros((nyears,12)) * np.nan
area_camscsm[0:85,:] = area_camscsm_init
volume_camscsm[0:85,:] = volume_camscsm_init

# Load sea-ice area and volume FGOALS-f3-L
filename = dir_input + 'SIA-SIV_SImon_FGOALS-f3-L_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_fgoalsf3l,volume_fgoalsf3l = np.load(filename)

# Load sea-ice area and volume FGOALS-g3
filename = dir_input + 'SIA-SIV_SImon_FGOALS-g3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_fgoalsg3,volume_fgoalsg3,notused,notused = np.load(filename)

# Load sea-ice area and volume CanESM5
filename = dir_input + 'SIA-SIV_SImon_CanESM5_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_canesm5,volume_canesm5,notused,notused = np.load(filename)

# Load sea-ice area and volume CanESM5-CanOE
filename = dir_input + 'SIA-SIV_SImon_CanESM5-CanOE_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_canesm5canoe,volume_canesm5canoe,notused,notused = np.load(filename)

# Load sea-ice area and volume CNRM-CM6-1
filename = dir_input + 'SIA-SIV_SImon_CNRM-CM6-1_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cnrmcm6,volume_cnrmcm6,notused,notused = np.load(filename)

# Load sea-ice area and volume CNRM-CM6-1-HR
filename = dir_input + 'SIA-SIV_SImon_CNRM-CM6-1-HR_' + str(experiment) + '_r1i1p1f2_gn_2015-2100.npy'
area_cnrmcm6hr,volume_cnrmcm6hr = np.load(filename)

# Load sea-ice area and volume ACCESS-ESM1-5
filename = dir_input + 'SIA-SIV_SImon_ACCESS-ESM1-5_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_accessesm,volume_accessesm,notused,notused = np.load(filename)

# Load sea-ice area and volume ACCESS-CM2
filename = dir_input + 'SIA-SIV_SImon_ACCESS-CM2_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_accesscm,volume_accesscm,notused,notused = np.load(filename)

# Load sea-ice area and volume MPI-ESM1-2-LR
filename = dir_input + 'SIA-SIV_SImon_MPI-ESM1-2-LR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_mpiesmlr,volume_mpiesmlr,notused,notused = np.load(filename)

# Load sea-ice area and volume MPI-ESM1-2-HR
filename = dir_input + 'SIA-SIV_SImon_MPI-ESM1-2-HR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_mpiesmhr,volume_mpiesmhr,notused,notused = np.load(filename)

# Load sea-ice area and volume EC-Earth3
filename = dir_input + 'SIA-SIV_SImon_EC-Earth3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ecearth3,volume_ecearth3,notused,notused = np.load(filename)

# Load sea-ice area and volume EC-Earth3-Veg
filename = dir_input + 'SIA-SIV_SImon_EC-Earth3-Veg_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ecearth3veg,volume_ecearth3veg,notused,notused = np.load(filename)

# Load sea-ice area and volume FIO-ESM-2-0
filename = dir_input + 'SIA-SIV_SImon_FIO-ESM-2-0_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_fioesm,volume_fioesm,notused,notused = np.load(filename)

# Load sea-ice area and volume INM-CM4-8
filename = dir_input + 'SIA-SIV_SImon_INM-CM4-8_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_inmcm48,volume_inmcm48 = np.load(filename)

# Load sea-ice area and volume INM-CM5-0
filename = dir_input + 'SIA-SIV_SImon_INM-CM5-0_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_inmcm50,volume_inmcm50 = np.load(filename)

# Load sea-ice area and volume IPSL-CM6A-LR
filename = dir_input + 'SIA-SIV_SImon_IPSL-CM6A-LR_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ipslcm6alr,volume_ipslcm6alr,notused,notused = np.load(filename)

# Load sea-ice area and volume KIOST-ESM
filename = dir_input + 'SIA-SIV_SImon_KIOST-ESM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_kiostesm,volume_kiostesm = np.load(filename)

# Load sea-ice area and volume MIROC6
filename = dir_input + 'SIA-SIV_SImon_MIROC6_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_miroc6,volume_miroc6,notused,notused = np.load(filename)

# Load sea-ice area and volume MIROC-ES2L
filename = dir_input + 'SIA-SIV_SImon_MIROC-ES2L_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_miroces2l,volume_miroces2l,notused,notused = np.load(filename)

# Load sea-ice area and volume HadGEM3-GC31-LL
filename = dir_input + 'SIA-SIV_SImon_HadGEM3-GC31-LL_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_hadgem3ll,volume_hadgem3ll,notused,notused = np.load(filename)

# Load sea-ice area and volume HadGEM3-GC31-MM
filename = dir_input + 'SIA-SIV_SImon_HadGEM3-GC31-MM_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_hadgem3mm,volume_hadgem3mm,notused,notused = np.load(filename)

# Load sea-ice area and volume UKESM1-0-LL
filename = dir_input + 'SIA-SIV_SImon_UKESM1-0-LL_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_ukesmll,volume_ukesmll,notused,notused = np.load(filename)

# Load sea-ice area and volume MRI-ESM2-0
filename = dir_input + 'SIA-SIV_SImon_MRI-ESM2-0_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_mriesm,volume_mriesm = np.load(filename)

# Load sea-ice area and volume CESM2
filename = dir_input + 'SIA-SIV_SImon_CESM2_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cesm2,volume_cesm2,notused,notused = np.load(filename)

# Load sea-ice area and volume CESM2-WACCM
filename = dir_input + 'SIA-SIV_SImon_CESM2-WACCM_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_cesm2waccm,volume_cesm2waccm,notused,notused = np.load(filename)

# Load sea-ice area and volume NorESM2-LM
filename = dir_input + 'SIA-SIV_SImon_NorESM2-LM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_noresm2lm,volume_noresm2lm = np.load(filename)

# Load sea-ice area and volume NorESM2-MM
filename = dir_input + 'SIA-SIV_SImon_NorESM2-MM_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_noresm2mm,volume_noresm2mm = np.load(filename)

# Load sea-ice area and volume GFDL-CM4
if experiment == 'ssp585':
    filename = dir_input + 'SIA-SIV_SImon_GFDL-CM4_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
    area_gfdlcm4,volume_gfdlcm4 = np.load(filename)

# Load sea-ice area and volume GFDL-ESM4
filename = dir_input + 'SIA-SIV_SImon_GFDL-ESM4_' + str(experiment) + '_r1i1p1f1_gn_2015-2100.npy'
area_gfdlesm4,volume_gfdlesm4 = np.load(filename)

# Load sea-ice area and volume NESM3
filename = dir_input + 'SIA-SIV_SImon_NESM3_' + str(experiment) + '_ensmean_gn_2015-2100.npy'
area_nesm3,volume_nesm3,notused,notused = np.load(filename)

# Compute multi-model mean - all models
area_mmm = np.zeros((nyears,nmy))
volume_mmm = np.zeros((nyears,nmy))
sd_area_mmm = np.zeros((nyears,nmy))
sd_volume_mmm = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        if experiment == 'ssp585':
            array_area = [area_awicm[y,m],area_bcccsm2mr[y,m],area_camscsm[y,m],area_fgoalsf3l[y,m],area_fgoalsg3[y,m],area_canesm5[y,m],area_canesm5canoe[y,m],area_cnrmcm6[y,m],area_cnrmcm6hr[y,m],area_accessesm[y,m],area_accesscm[y,m],area_mpiesmlr[y,m],area_mpiesmhr[y,m],area_ecearth3[y,m],area_ecearth3veg[y,m],area_fioesm[y,m],area_inmcm48[y,m],area_inmcm50[y,m],area_ipslcm6alr[y,m],area_kiostesm[y,m],area_miroc6[y,m],area_miroces2l[y,m],area_hadgem3ll[y,m],area_hadgem3mm[y,m],area_ukesmll[y,m],area_mriesm[y,m],area_cesm2[y,m],area_cesm2waccm[y,m],area_noresm2lm[y,m],area_noresm2mm[y,m],area_gfdlcm4[y,m],area_gfdlesm4[y,m],area_nesm3[y,m]]
            array_volume = [volume_awicm[y,m],volume_bcccsm2mr[y,m],volume_camscsm[y,m],volume_fgoalsf3l[y,m],volume_canesm5[y,m],volume_cnrmcm6[y,m],volume_cnrmcm6hr[y,m],volume_accessesm[y,m],volume_accesscm[y,m],volume_mpiesmlr[y,m],volume_mpiesmhr[y,m],volume_ecearth3[y,m],volume_ecearth3veg[y,m],volume_fioesm[y,m],volume_ipslcm6alr[y,m],volume_kiostesm[y,m],volume_miroc6[y,m],volume_miroces2l[y,m],volume_hadgem3ll[y,m],volume_hadgem3mm[y,m],volume_ukesmll[y,m],volume_mriesm[y,m],volume_cesm2[y,m],volume_cesm2waccm[y,m],volume_noresm2lm[y,m],volume_noresm2mm[y,m],volume_gfdlcm4[y,m],volume_nesm3[y,m]]
        else:
            array_area = [area_awicm[y,m],area_bcccsm2mr[y,m],area_camscsm[y,m],area_fgoalsf3l[y,m],area_fgoalsg3[y,m],area_canesm5[y,m],area_canesm5canoe[y,m],area_cnrmcm6[y,m],area_cnrmcm6hr[y,m],area_accessesm[y,m],area_accesscm[y,m],area_mpiesmlr[y,m],area_mpiesmhr[y,m],area_ecearth3[y,m],area_ecearth3veg[y,m],area_fioesm[y,m],area_inmcm48[y,m],area_inmcm50[y,m],area_ipslcm6alr[y,m],area_kiostesm[y,m],area_miroc6[y,m],area_miroces2l[y,m],area_hadgem3ll[y,m],area_hadgem3mm[y,m],area_ukesmll[y,m],area_mriesm[y,m],area_cesm2[y,m],area_cesm2waccm[y,m],area_noresm2lm[y,m],area_noresm2mm[y,m],area_gfdlesm4[y,m],area_nesm3[y,m]]
            array_volume = [volume_awicm[y,m],volume_bcccsm2mr[y,m],volume_camscsm[y,m],volume_fgoalsf3l[y,m],volume_canesm5[y,m],volume_cnrmcm6[y,m],volume_cnrmcm6hr[y,m],volume_accessesm[y,m],volume_accesscm[y,m],volume_mpiesmlr[y,m],volume_mpiesmhr[y,m],volume_ecearth3[y,m],volume_ecearth3veg[y,m],volume_fioesm[y,m],volume_ipslcm6alr[y,m],volume_kiostesm[y,m],volume_miroc6[y,m],volume_miroces2l[y,m],volume_hadgem3ll[y,m],volume_hadgem3mm[y,m],volume_ukesmll[y,m],volume_mriesm[y,m],volume_cesm2[y,m],volume_cesm2waccm[y,m],volume_noresm2lm[y,m],volume_noresm2mm[y,m],volume_gfdlesm4[y,m],volume_nesm3[y,m]]
        area_mmm[y,m] = np.nanmean(array_area)
        volume_mmm[y,m] = np.nanmean(array_volume)
        sd_area_mmm[y,m] = np.nanstd(array_area)
        sd_volume_mmm[y,m] = np.nanstd(array_volume)

# Compute multi-model mean - good mean SIA (15 best models)
area_select_sia_15 = np.zeros((nyears,nmy))
volume_select_sia_15 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        if experiment == 'ssp585':
            array_area = [area_cesm2waccm[y,m],area_noresm2lm[y,m],area_accessesm[y,m],area_gfdlesm4[y,m],area_mpiesmlr[y,m],area_cnrmcm6hr[y,m],area_hadgem3mm[y,m],area_ipslcm6alr[y,m],area_gfdlcm4[y,m],area_inmcm50[y,m],area_cnrmcm6[y,m],area_ecearth3[y,m],area_ecearth3veg[y,m],area_hadgem3ll[y,m],area_miroces2l[y,m]]
            array_volume = [volume_cesm2waccm[y,m],volume_noresm2lm[y,m],volume_accessesm[y,m],volume_gfdlesm4[y,m],volume_mpiesmlr[y,m],volume_cnrmcm6hr[y,m],volume_hadgem3mm[y,m],volume_ipslcm6alr[y,m],volume_gfdlcm4[y,m],volume_cnrmcm6[y,m],volume_ecearth3[y,m],volume_ecearth3veg[y,m],volume_hadgem3ll[y,m],volume_miroces2l[y,m],volume_mriesm[y,m]]
        else:
            array_area = [area_cesm2waccm[y,m],area_noresm2lm[y,m],area_accessesm[y,m],area_gfdlesm4[y,m],area_mpiesmlr[y,m],area_cnrmcm6hr[y,m],area_hadgem3mm[y,m],area_ipslcm6alr[y,m],area_inmcm50[y,m],area_cnrmcm6[y,m],area_ecearth3[y,m],area_ecearth3veg[y,m],area_hadgem3ll[y,m],area_miroces2l[y,m],area_mriesm[y,m]]
            array_volume = [volume_cesm2waccm[y,m],volume_noresm2lm[y,m],volume_accessesm[y,m],volume_gfdlesm4[y,m],volume_mpiesmlr[y,m],volume_cnrmcm6hr[y,m],volume_hadgem3mm[y,m],volume_ipslcm6alr[y,m],volume_cnrmcm6[y,m],volume_ecearth3[y,m],volume_ecearth3veg[y,m],volume_hadgem3ll[y,m],volume_miroces2l[y,m],volume_mriesm[y,m],volume_accesscm[y,m]]
        area_select_sia_15[y,m] = np.nanmean(array_area)
        volume_select_sia_15[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good mean SIV (15 best models)
area_select_siv_15 = np.zeros((nyears,nmy))
volume_select_siv_15 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        if experiment == 'ssp585':
            array_area = [area_cesm2waccm[y,m],area_ipslcm6alr[y,m],area_fgoalsf3l[y,m],area_camscsm[y,m],area_miroc6[y,m],area_mpiesmhr[y,m],area_fioesm[y,m],area_accessesm[y,m],area_mpiesmlr[y,m],area_awicm[y,m],area_cesm2[y,m],area_gfdlcm4[y,m],area_hadgem3mm[y,m],area_mriesm[y,m],area_accesscm[y,m]]
            array_volume = [volume_cesm2waccm[y,m],volume_ipslcm6alr[y,m],volume_fgoalsf3l[y,m],volume_camscsm[y,m],volume_miroc6[y,m],volume_mpiesmhr[y,m],volume_fioesm[y,m],volume_accessesm[y,m],volume_mpiesmlr[y,m],volume_awicm[y,m],volume_cesm2[y,m],volume_gfdlcm4[y,m],volume_hadgem3mm[y,m],volume_mriesm[y,m],volume_accesscm[y,m]]
        else:
            array_area = [area_cesm2waccm[y,m],area_ipslcm6alr[y,m],area_fgoalsf3l[y,m],area_camscsm[y,m],area_miroc6[y,m],area_mpiesmhr[y,m],area_fioesm[y,m],area_accessesm[y,m],area_mpiesmlr[y,m],area_awicm[y,m],area_cesm2[y,m],area_hadgem3mm[y,m],area_mriesm[y,m],area_accesscm[y,m],area_noresm2lm[y,m]]
            array_volume = [volume_cesm2waccm[y,m],volume_ipslcm6alr[y,m],volume_fgoalsf3l[y,m],volume_camscsm[y,m],volume_miroc6[y,m],volume_mpiesmhr[y,m],volume_fioesm[y,m],volume_accessesm[y,m],volume_mpiesmlr[y,m],volume_awicm[y,m],volume_cesm2[y,m],volume_hadgem3mm[y,m],volume_mriesm[y,m],volume_accesscm[y,m],volume_noresm2lm[y,m]]
        area_select_siv_15[y,m] = np.nanmean(array_area)
        volume_select_siv_15[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good trend in SIA (15 best models)
area_select_tsia_15 = np.zeros((nyears,nmy))
volume_select_tsia_15 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        if experiment == 'ssp585':
            array_area = [area_hadgem3mm[y,m],area_cesm2waccm[y,m],area_fgoalsf3l[y,m],area_mpiesmhr[y,m],area_cesm2[y,m],area_mriesm[y,m],area_accessesm[y,m],area_bcccsm2mr[y,m],area_fioesm[y,m],area_gfdlesm4[y,m],area_ipslcm6alr[y,m],area_mpiesmlr[y,m],area_cnrmcm6hr[y,m],area_miroc6[y,m],area_gfdlcm4[y,m]]
            array_volume = [volume_hadgem3mm[y,m],volume_cesm2waccm[y,m],volume_fgoalsf3l[y,m],volume_mpiesmhr[y,m],volume_cesm2[y,m],volume_mriesm[y,m],volume_accessesm[y,m],volume_bcccsm2mr[y,m],volume_fioesm[y,m],volume_gfdlesm4[y,m],volume_ipslcm6alr[y,m],volume_mpiesmlr[y,m],volume_cnrmcm6hr[y,m],volume_miroc6[y,m],volume_gfdlcm4[y,m]]
        else:
            array_area = [area_hadgem3mm[y,m],area_cesm2waccm[y,m],area_fgoalsf3l[y,m],area_mpiesmhr[y,m],area_cesm2[y,m],area_mriesm[y,m],area_accessesm[y,m],area_bcccsm2mr[y,m],area_fioesm[y,m],area_gfdlesm4[y,m],area_ipslcm6alr[y,m],area_mpiesmlr[y,m],area_cnrmcm6hr[y,m],area_miroc6[y,m],area_hadgem3ll[y,m]]
            array_volume = [volume_hadgem3mm[y,m],volume_cesm2waccm[y,m],volume_fgoalsf3l[y,m],volume_mpiesmhr[y,m],volume_cesm2[y,m],volume_mriesm[y,m],volume_accessesm[y,m],volume_bcccsm2mr[y,m],volume_fioesm[y,m],volume_gfdlesm4[y,m],volume_ipslcm6alr[y,m],volume_mpiesmlr[y,m],volume_cnrmcm6hr[y,m],volume_miroc6[y,m],volume_hadgem3ll[y,m]]
        area_select_tsia_15[y,m] = np.nanmean(array_area)
        volume_select_tsia_15[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good trend in SIV (15 best models)
area_select_tsiv_15 = np.zeros((nyears,nmy))
volume_select_tsiv_15 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        if experiment == 'ssp585':
            array_area = [area_ipslcm6alr[y,m],area_cesm2[y,m],area_fioesm[y,m],area_nesm3[y,m],area_noresm2mm[y,m],area_noresm2lm[y,m],area_mpiesmhr[y,m],area_miroc6[y,m],area_accesscm[y,m],area_gfdlcm4[y,m],area_mpiesmlr[y,m],area_fgoalsf3l[y,m],area_gfdlesm4[y,m],area_kiostesm[y,m],area_ecearth3[y,m]]
            array_volume = [volume_ipslcm6alr[y,m],volume_cesm2[y,m],volume_fioesm[y,m],volume_nesm3[y,m],volume_noresm2mm[y,m],volume_noresm2lm[y,m],volume_mpiesmhr[y,m],volume_miroc6[y,m],volume_accesscm[y,m],volume_gfdlcm4[y,m],volume_mpiesmlr[y,m],volume_fgoalsf3l[y,m],volume_gfdlesm4[y,m],volume_kiostesm[y,m],volume_ecearth3[y,m]]
        else:
            array_area = [area_ipslcm6alr[y,m],area_cesm2[y,m],area_fioesm[y,m],area_nesm3[y,m],area_noresm2mm[y,m],area_noresm2lm[y,m],area_mpiesmhr[y,m],area_miroc6[y,m],area_accesscm[y,m],area_mpiesmlr[y,m],area_fgoalsf3l[y,m],area_gfdlesm4[y,m],area_kiostesm[y,m],area_ecearth3[y,m],area_mriesm[y,m]]
            array_volume = [volume_ipslcm6alr[y,m],volume_cesm2[y,m],volume_fioesm[y,m],volume_nesm3[y,m],volume_noresm2mm[y,m],volume_noresm2lm[y,m],volume_mpiesmhr[y,m],volume_miroc6[y,m],volume_accesscm[y,m],volume_mpiesmlr[y,m],volume_fgoalsf3l[y,m],volume_gfdlesm4[y,m],volume_kiostesm[y,m],volume_ecearth3[y,m],volume_mriesm[y,m]]
        area_select_tsiv_15[y,m] = np.nanmean(array_area)
        volume_select_tsiv_15[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good SIA variability (15 best models)
area_select_varsia_15 = np.zeros((nyears,nmy))
volume_select_varsia_15 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        if experiment == 'ssp585':
            array_area = [area_gfdlcm4[y,m],area_accessesm[y,m],area_inmcm50[y,m],area_noresm2mm[y,m],area_gfdlesm4[y,m],area_canesm5canoe[y,m],area_cnrmcm6hr[y,m],area_accesscm[y,m],area_hadgem3ll[y,m],area_nesm3[y,m],area_fgoalsf3l[y,m],area_fgoalsg3[y,m],area_inmcm48[y,m],area_bcccsm2mr[y,m],area_ecearth3[y,m]]
            array_volume = [volume_gfdlcm4[y,m],volume_accessesm[y,m],volume_noresm2mm[y,m],volume_gfdlesm4[y,m],volume_cnrmcm6hr[y,m],volume_accesscm[y,m],volume_hadgem3ll[y,m],volume_nesm3[y,m],volume_fgoalsf3l[y,m],volume_bcccsm2mr[y,m],volume_ecearth3[y,m],volume_ecearth3veg[y,m],volume_cesm2waccm[y,m],volume_kiostesm[y,m],volume_camscsm[y,m]]
        else:
            array_area = [area_accessesm[y,m],area_inmcm50[y,m],area_noresm2mm[y,m],area_gfdlesm4[y,m],area_canesm5canoe[y,m],area_cnrmcm6hr[y,m],area_accesscm[y,m],area_hadgem3ll[y,m],area_nesm3[y,m],area_fgoalsf3l[y,m],area_fgoalsg3[y,m],area_inmcm48[y,m],area_bcccsm2mr[y,m],area_ecearth3[y,m],area_mriesm[y,m]]
            array_volume = [volume_accessesm[y,m],volume_noresm2mm[y,m],volume_gfdlesm4[y,m],volume_cnrmcm6hr[y,m],volume_accesscm[y,m],volume_hadgem3ll[y,m],volume_nesm3[y,m],volume_fgoalsf3l[y,m],volume_bcccsm2mr[y,m],volume_ecearth3[y,m],volume_ecearth3veg[y,m],volume_cesm2waccm[y,m],volume_kiostesm[y,m],volume_camscsm[y,m],volume_mriesm[y,m]] 
        area_select_varsia_15[y,m] = np.nanmean(array_area)
        volume_select_varsia_15[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good SIV variability (15 best models)
area_select_varsiv_15 = np.zeros((nyears,nmy))
volume_select_varsiv_15 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        if experiment == 'ssp585':
            array_area = [area_gfdlcm4[y,m],area_canesm5[y,m],area_noresm2lm[y,m],area_noresm2mm[y,m],area_fioesm[y,m],area_hadgem3mm[y,m],area_cesm2waccm[y,m],area_miroces2l[y,m],area_accesscm[y,m],area_gfdlesm4[y,m],area_cesm2[y,m],area_ecearth3veg[y,m],area_ecearth3[y,m],area_fgoalsf3l[y,m],area_camscsm[y,m]]
            array_volume = [volume_gfdlcm4[y,m],volume_canesm5[y,m],volume_noresm2lm[y,m],volume_noresm2mm[y,m],volume_fioesm[y,m],volume_hadgem3mm[y,m],volume_cesm2waccm[y,m],volume_miroces2l[y,m],volume_accesscm[y,m],volume_gfdlesm4[y,m],volume_cesm2[y,m],volume_ecearth3veg[y,m],volume_ecearth3[y,m],volume_fgoalsf3l[y,m],volume_camscsm[y,m]]
        else:
            array_area = [area_canesm5[y,m],area_noresm2lm[y,m],area_noresm2mm[y,m],area_fioesm[y,m],area_hadgem3mm[y,m],area_cesm2waccm[y,m],area_miroces2l[y,m],area_accesscm[y,m],area_gfdlesm4[y,m],area_cesm2[y,m],area_ecearth3veg[y,m],area_ecearth3[y,m],area_fgoalsf3l[y,m],area_camscsm[y,m],area_ukesmll[y,m]]
            array_volume = [volume_canesm5[y,m],volume_noresm2lm[y,m],volume_noresm2mm[y,m],volume_fioesm[y,m],volume_hadgem3mm[y,m],volume_cesm2waccm[y,m],volume_miroces2l[y,m],volume_accesscm[y,m],volume_gfdlesm4[y,m],volume_cesm2[y,m],volume_ecearth3veg[y,m],volume_ecearth3[y,m],volume_fgoalsf3l[y,m],volume_camscsm[y,m],volume_ukesmll[y,m]]
        area_select_varsiv_15[y,m] = np.nanmean(array_area)
        volume_select_varsiv_15[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good AOHT 26N and 57N (8 best models)
area_select_aoht_8 = np.zeros((nyears,nmy))
volume_select_aoht_8 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        array_area = [area_hadgem3mm[y,m],area_ecearth3veg[y,m],area_mriesm[y,m],area_mpiesmlr[y,m],area_ecearth3[y,m],area_hadgem3ll[y,m],area_ukesmll[y,m],area_cnrmcm6[y,m]]
        array_volume = [volume_hadgem3mm[y,m],volume_ecearth3veg[y,m],volume_mriesm[y,m],volume_mpiesmlr[y,m],volume_ecearth3[y,m],volume_hadgem3ll[y,m],volume_ukesmll[y,m],volume_cnrmcm6[y,m]]
        area_select_aoht_8[y,m] = np.nanmean(array_area)
        volume_select_aoht_8[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good AOHT 70N and POHT 60N (8 best models)
area_select_apoht_8 = np.zeros((nyears,nmy))
volume_select_apoht_8 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        array_area = [area_mpiesmlr[y,m],area_ukesmll[y,m],area_mpiesmhr[y,m],area_hadgem3ll[y,m],area_hadgem3mm[y,m],area_ipslcm6alr[y,m],area_mriesm[y,m],area_canesm5[y,m]]
        array_volume = [volume_mpiesmlr[y,m],volume_ukesmll[y,m],volume_mpiesmhr[y,m],volume_hadgem3ll[y,m],volume_hadgem3mm[y,m],volume_ipslcm6alr[y,m],volume_mriesm[y,m],volume_canesm5[y,m]]
        area_select_apoht_8[y,m] = np.nanmean(array_area)
        volume_select_apoht_8[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good AOHT 26N+57N and SIA (among the best 15 SIA models and best 8 OHT models)
area_select_ohtsia1 = np.zeros((nyears,nmy))
volume_select_ohtsia1 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        array_area = [area_hadgem3mm[y,m],area_ecearth3veg[y,m],area_mpiesmlr[y,m],area_ecearth3[y,m],area_hadgem3ll[y,m]]
        array_volume = [volume_hadgem3mm[y,m],volume_ecearth3veg[y,m],volume_mpiesmlr[y,m],volume_ecearth3[y,m],volume_hadgem3ll[y,m]]
        area_select_ohtsia1[y,m] = np.nanmean(array_area)
        volume_select_ohtsia1[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good AOHT 70N + POHT 60N and SIA (among the best 15 SIA models and best 8 OHT models)
area_select_ohtsia2 = np.zeros((nyears,nmy))
volume_select_ohtsia2 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        array_area = [area_mpiesmlr[y,m],area_hadgem3ll[y,m],area_hadgem3mm[y,m],area_ipslcm6alr[y,m],area_canesm5[y,m]]
        array_volume = [volume_mpiesmlr[y,m],volume_hadgem3ll[y,m],volume_hadgem3mm[y,m],volume_ipslcm6alr[y,m],volume_canesm5[y,m]]
        area_select_ohtsia2[y,m] = np.nanmean(array_area)
        volume_select_ohtsia2[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good AOHT 26N+57N and SIV (among the best 15 SIV models and best 8 OHT models)
area_select_ohtsiv1 = np.zeros((nyears,nmy))
volume_select_ohtsiv1 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        array_area = [area_hadgem3mm[y,m],area_mriesm[y,m],area_mpiesmlr[y,m]]
        array_volume = [volume_hadgem3mm[y,m],volume_mriesm[y,m],volume_mpiesmlr[y,m]]
        area_select_ohtsiv1[y,m] = np.nanmean(array_area)
        volume_select_ohtsiv1[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - good AOHT 70N + POHT 60N and SIV (among the best 15 SIV models and best 8 OHT models)
area_select_ohtsiv2 = np.zeros((nyears,nmy))
volume_select_ohtsiv2 = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        array_area = [area_mpiesmlr[y,m],area_mpiesmhr[y,m],area_hadgem3mm[y,m],area_ipslcm6alr[y,m],area_mriesm[y,m],area_canesm5[y,m]]
        array_volume = [volume_mpiesmlr[y,m],volume_mpiesmhr[y,m],volume_hadgem3mm[y,m],volume_ipslcm6alr[y,m],volume_mriesm[y,m],volume_canesm5[y,m]]
        area_select_ohtsiv2[y,m] = np.nanmean(array_area)
        volume_select_ohtsiv2[y,m] = np.nanmean(array_volume)

# Compute multi-model mean - >=5 members
area_select_members = np.zeros((nyears,nmy))
volume_select_members = np.zeros((nyears,nmy))
for m in np.arange(nmy):
    for y in np.arange(nyears):
        if experiment == 'ssp585':
            array_area = [area_canesm5[y,m],area_cnrmcm6[y,m],area_mpiesmlr[y,m],area_ecearth3[y,m],area_ecearth3veg[y,m],area_ipslcm6alr[y,m],area_miroc6[y,m],area_ukesmll[y,m],area_cesm2[y,m],area_cesm2waccm[y,m]]
            array_volume = [volume_canesm5[y,m],volume_cnrmcm6[y,m],volume_mpiesmlr[y,m],volume_ecearth3[y,m],volume_ecearth3veg[y,m],volume_ipslcm6alr[y,m],volume_miroc6[y,m],volume_ukesmll[y,m],volume_cesm2[y,m],volume_cesm2waccm[y,m]]
        elif experiment == 'ssp126':
            array_area = [area_canesm5[y,m],area_cnrmcm6[y,m],area_mpiesmlr[y,m],area_ecearth3[y,m],area_ecearth3veg[y,m],area_ipslcm6alr[y,m],area_miroc6[y,m],area_ukesmll[y,m],area_cesm2[y,m]]
            array_volume = [volume_canesm5[y,m],volume_cnrmcm6[y,m],volume_mpiesmlr[y,m],volume_ecearth3[y,m],volume_ecearth3veg[y,m],volume_ipslcm6alr[y,m],volume_miroc6[y,m],volume_ukesmll[y,m],volume_cesm2[y,m]]
        area_select_members[y,m] = np.nanmean(array_area)
        volume_select_members[y,m] = np.nanmean(array_volume)


# Print
print(np.arange(nyears)+2015)
print(area_select_ohtsiv2[:,8])

# Labels
name_xticks = ['2010','2020','2030','2040','2050','2060','2070','2080','2090','2100']
    
# Time series of total Arctic sea-ice area
xrangea = np.arange(1, 92, 10)
fig,ax = plt.subplots(2,1,figsize=(20,12))
fig.subplots_adjust(left=0.08,bottom=0.09,right=0.7,top=0.95,wspace=None,hspace=0.3)

# March area
if experiment == 'ssp126':
    ax[0].plot(np.arange(nyears) + 6,area_mmm[:,2],'-',color='black',label='Without selection (32)',linewidth=4)
elif experiment == 'ssp585':
    ax[0].plot(np.arange(nyears) + 6,area_mmm[:,2],'-',color='black',label='Without selection (33)',linewidth=4)
ax[0].fill_between(np.arange(nyears) + 6,area_mmm[:,2]-sd_area_mmm[:,2],area_mmm[:,2]+sd_area_mmm[:,2],alpha=0.5,edgecolor='gray',facecolor='gray',linestyle='--')
ax[0].plot(np.arange(nyears) + 6,area_select_sia_15[:,2],'-',color='blue',label='Mean sea-ice area (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_siv_15[:,2],'-',color='green',label='Mean sea-ice volume (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_varsia_15[:,2],'-.',color='blue',label='Sea-ice area variability (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_varsiv_15[:,2],'-.',color='green',label='Sea-ice volume variability (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_tsia_15[:,2],'--',color='blue',label='Trend in sea-ice area (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_tsiv_15[:,2],'--',color='green',label='Trend in sea-ice volume (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_aoht_8[:,2],'-',color='red',label='Atlantic OHT (8)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_apoht_8[:,2],'--',color='red',label='Atlantic/Pacific OHT (8)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_ohtsia1[:,2],'-',color='orange',label='Atlantic OHT + sea-ice area (5)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_ohtsia2[:,2],'--',color='orange',label='Atl/Pac OHT + sea-ice area (5)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_ohtsiv1[:,2],'-',color='gray',label='Atlantic OHT + sea-ice volume (3)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,area_select_ohtsiv2[:,2],'--',color='gray',label='Atl/Pac OHT + sea-ice volume (6)',linewidth=2)
if experiment == 'ssp585':
    ax[0].plot(np.arange(nyears) + 6,area_select_members[:,2],'-',color='magenta',label='>= 5 members (10)',linewidth=2)
elif experiment == 'ssp126':
    ax[0].plot(np.arange(nyears) + 6,area_select_members[:,2],'-',color='magenta',label='>= 5 members (9)',linewidth=2)
ax[0].legend(shadow=True,frameon=False,fontsize=20,bbox_to_anchor=(1,1))
ax[0].set_ylabel('March sea-ice area (10$^6$ km$^2$)',fontsize=24)
ax[0].set_xlabel('Year',fontsize=24)
ax[0].set_xticks(xrangea)
ax[0].set_xticklabels(name_xticks)
ax[0].set_yticks(np.arange(2, 16.1, 2))
ax[0].tick_params(axis='both',labelsize=20)
ax[0].axis([-1, 93, 2, 16])
ax[0].grid(linestyle='--')
ax[0].set_title('a',loc='left',fontsize=25,fontweight='bold')
    
# September area
ax[1].plot(np.arange(nyears) + 6,area_mmm[:,8],'-',color='black',linewidth=4)
ax[1].fill_between(np.arange(nyears) + 6,area_mmm[:,8]-sd_area_mmm[:,8],area_mmm[:,8]+sd_area_mmm[:,8],alpha=0.5,edgecolor='gray',facecolor='gray',linestyle='--')
ax[1].plot(np.arange(nyears) + 6,area_select_sia_15[:,8],'-',color='blue',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_siv_15[:,8],'-',color='green',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_varsia_15[:,8],'-.',color='blue',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_varsiv_15[:,8],'-.',color='green',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_tsia_15[:,8],'--',color='blue',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_tsiv_15[:,8],'--',color='green',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_aoht_8[:,8],'-',color='red',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_apoht_8[:,8],'--',color='red',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_ohtsia1[:,8],'-',color='orange',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_ohtsia2[:,8],'--',color='orange',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_ohtsiv1[:,8],'-',color='gray',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_ohtsiv2[:,8],'--',color='gray',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,area_select_members[:,8],'-',color='magenta',linewidth=2)
plt.axhline(y=1,color='black',linestyle='--',linewidth=2)
ax[1].set_xlabel('Year',fontsize=24)
ax[1].set_ylabel('September sea-ice area (10$^6$ km$^2$)',fontsize=24)
ax[1].set_xticks(xrangea)
ax[1].set_xticklabels(name_xticks)
ax[1].set_yticks(np.arange(0, 6.1, 1))
ax[1].tick_params(axis='both',labelsize=20)
ax[1].axis([-1, 93, 0, 6])
ax[1].grid(linestyle='--')
ax[1].set_title('b',loc='left',fontsize=25,fontweight='bold')
    
if save_fig == True:
    if experiment == 'ssp585':
        fig.savefig(dir_output + 'fig2.png')
    else:
        fig.savefig(dir_output + 'figS2.png')
    
    
# Time series of total Arctic sea-ice volume
xrangea = np.arange(1, 92, 10)
fig,ax = plt.subplots(2,1,figsize=(20,12))
fig.subplots_adjust(left=0.08,bottom=0.09,right=0.7,top=0.95,wspace=None,hspace=0.3)
   
# March volume
ax[0].plot(np.arange(nyears) + 6,volume_mmm[:,2],'-',color='black',label='Without selection (28)',linewidth=3)
ax[0].fill_between(np.arange(nyears) + 6,volume_mmm[:,2]-sd_volume_mmm[:,2],volume_mmm[:,2]+sd_volume_mmm[:,2],alpha=0.5,edgecolor='gray',facecolor='gray',linestyle='--')
ax[0].plot(np.arange(nyears) + 6,volume_select_sia_15[:,2],'-',color='blue',label='Mean sea-ice area (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_siv_15[:,2],'-',color='green',label='Mean sea-ice volume (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_varsia_15[:,2],'-.',color='blue',label='Sea-ice area variability (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_varsiv_15[:,2],'-.',color='green',label='Sea-ice volume variability (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_tsia_15[:,2],'--',color='blue',label='Trend in sea-ice area (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_tsiv_15[:,2],'--',color='green',label='Trend in sea-ice volume (15)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_aoht_8[:,2],'-',color='red',label='Atlantic OHT (8)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_apoht_8[:,2],'--',color='red',label='Atlantic/Pacific OHT (8)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_ohtsia1[:,2],'-',color='orange',label='Atlantic OHT and SIA (5)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_ohtsia2[:,2],'--',color='orange',label='Atl/Pac OHT and sea-ice area (5)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_ohtsiv1[:,2],'-',color='gray',label='Atlantic OHT and sea-ice volume (3)',linewidth=2)
ax[0].plot(np.arange(nyears) + 6,volume_select_ohtsiv2[:,2],'--',color='gray',label='Atl/Pac OHT and sea-ice volume (6)',linewidth=2)
if experiment == 'ssp585':
    ax[0].plot(np.arange(nyears) + 6,volume_select_members[:,2],'-',color='magenta',label='>= 5 members (10)',linewidth=2)
elif experiment == 'ssp126':
    ax[0].plot(np.arange(nyears) + 6,volume_select_members[:,2],'-',color='magenta',label='>= 5 members (9)',linewidth=2)
ax[0].legend(shadow=True,frameon=False,fontsize=20,bbox_to_anchor=(1,1))
ax[0].set_ylabel('March SIV (10$^3$ km$^3$)',fontsize=24)
ax[0].set_xticks(xrangea)
ax[0].set_xticklabels(name_xticks)
ax[0].set_yticks(np.arange(0, 30.1, 5))
ax[0].tick_params(axis='both',labelsize=20)
ax[0].axis([-1, 93, 0, 30])
ax[0].grid(linestyle='--')
ax[0].set_title('a',loc='left',fontsize=25,fontweight='bold')
    
# September volume
ax[1].plot(np.arange(nyears) + 6,volume_mmm[:,8],'-',color='black',linewidth=3)
ax[1].fill_between(np.arange(nyears) + 6,volume_mmm[:,8]-sd_volume_mmm[:,8],volume_mmm[:,8]+sd_volume_mmm[:,8],alpha=0.5,edgecolor='gray',facecolor='gray',linestyle='--')
ax[1].plot(np.arange(nyears) + 6,volume_select_sia_15[:,8],'-',color='blue',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_siv_15[:,8],'-',color='green',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_varsia_15[:,8],'-.',color='blue',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_varsiv_15[:,8],'-.',color='green',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_tsia_15[:,8],'--',color='blue',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_tsiv_15[:,8],'--',color='green',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_aoht_8[:,8],'-',color='red',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_apoht_8[:,8],'--',color='red',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_ohtsia1[:,8],'-',color='orange',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_ohtsia2[:,8],'--',color='orange',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_ohtsiv1[:,8],'-',color='gray',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_ohtsiv2[:,8],'--',color='gray',linewidth=2)
ax[1].plot(np.arange(nyears) + 6,volume_select_members[:,8],'-',color='magenta',linewidth=2)
ax[1].set_xlabel('Year',fontsize=24)
ax[1].set_ylabel('September SIV (10$^3$ km$^3$)',fontsize=24)
ax[1].set_xticks(xrangea)
ax[1].set_xticklabels(name_xticks)
ax[1].set_yticks(np.arange(0, 12.1, 2))
ax[1].tick_params(axis='both',labelsize=20)
ax[1].axis([-1, 93, 0, 12])
ax[1].grid(linestyle='--')
ax[1].set_title('b',loc='left',fontsize=25,fontweight='bold')

if save_fig == True:
    if experiment == 'ssp585':
        fig.savefig(dir_output + 'figS3.png')
    else:
        fig.savefig(dir_output + 'figS4.png')
