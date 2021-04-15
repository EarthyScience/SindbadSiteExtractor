#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 14:49:25 2019

@author: simon
"""
#%% Load library
import xarray as xr
import matplotlib.pyplot as plt
import os
from datetime import datetime
import json
import argparse
import sys
import numpy as np
sys.path.append('/home/simon/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/fluxnet_workflow/utils')
from utils.data_providers import EddyCovariance_FLUXNET, SWC_TS_FLUXNET, Climate_FLUXNET, Climate_ERA5, gapfill_climate, EvapTrans_Nelson2019, MODIS_MCD43A, MODIS_MxD11A, insitu_AGB, gloBbiomass_AGB, hansen_treecover, soilGrids_soilTexture, forest_age, insitu_VegFrac, hilda_VegFrac, data_harm
CurrentScript = os.path.basename('/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/fluxnet_workflow/fluxnet_workflow.py')

#%% Retrieve arguments to parse in data processing routine
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config_file',
                        '-config_file',
                        type=str,
                        help='Path to the config file',
                        default='/home/simon/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/fluxnet_workflow/exp_config.json',
                        required=False)   
    args = parser.parse_args()

    return args

#%% Retrieve parsed arguments
args = parse_args()

#%% Open config files
with open(args.config_file, 'r') as myfile:
    exp_config=myfile.read()
exp_config = json.loads(exp_config)

#%% Check sites to process
if exp_config["sites_toRun"] == 'all':
    sites = xr.open_zarr(exp_config["DataPath"]['fluxcom_cube']).site.values
else:
    sites = exp_config["sites_toRun"]

#%% Run data processing for each site
for site in sites:
    
    print(f'Processing data for site {site}')    
    
    #%% Process eddy-covariance data
    EddyCovariance_FLUXNET_data = EddyCovariance_FLUXNET(site= site, config=exp_config).EddyCovariance_FLUXNET_proc()
    
    #%% Process soil water content and soil temperature data
    SWC_TS_data = SWC_TS_FLUXNET(site= site, config=exp_config).SWC_TS_FLUXNET_proc()    
    
    #%% Process FLUXNET climate data
    Climate_FLUXNET_data = Climate_FLUXNET(site= site, config=exp_config).Climate_FLUXNET_proc()
        
    #%% Process downscaled climate data - ERA5
    Climate_ERA5_data = Climate_ERA5(site = site, config=exp_config).Climate_ERA5_proc()
            
    #%% Create gapfilled climate data
    gapfill_climate_data = gapfill_climate(Climate_site_data = Climate_FLUXNET_data, Climate_downscaled_data = Climate_ERA5_data, config = exp_config).gapfill_climate_proc()
    
    #%% Add T and E Nelson et al. data
    TEA_T_E = EvapTrans_Nelson2019(site = site, config=exp_config).EvapTrans_Nelson2019_proc()    
    
    #%% Process MODIS MCD43A data
    MCD43A_data = MODIS_MCD43A(site = site, config=exp_config).MODIS_MCD43A_proc()    
    
    #%% Process MODIS MxD11A data
    MxD11A_data = MODIS_MxD11A(site = site, config=exp_config).MODIS_MxD11A_proc()        
    
    #%% Process in-situ AGB data
    insitu_AGB_data = insitu_AGB(site = site, config=exp_config).insitu_AGB_proc()
    
    #%% Process globbiomass AGB data
    gloBbiomass_AGB_data = gloBbiomass_AGB(site = site, config=exp_config).gloBbiomass_AGB_proc()
    
    #%% Process Hansen tree cover data
    hansen_treecover_data = hansen_treecover(site = site, config=exp_config).hansen_treecover_proc() 
    
    #%% Process Soilgrids data
    soilGrids_soilTexture_data = soilGrids_soilTexture(site = site, config=exp_config).soilGrids_soilTexture_proc() 
        
    #%% Process age data
    DIST_FOREST_FRAC, LAST_DISTURBANCE_DATE = forest_age(site = site, config=exp_config).forest_age_proc()
    
    #%% Process insitu PFT variables
    PFT_, Frac_tree_insitu, Frac_grass_insitu, Frac_shrub_insitu, Frac_crop_insitu, Frac_savanna_insitu, Frac_wetland_insitu = insitu_VegFrac(site = site, config=exp_config).insitu_VegFrac_proc()    
    
    #%% Process hilda+data
    Frac_tree_hilda, Frac_grass_hilda, Frac_shrub_hilda, Frac_crop_hilda, Frac_other_hilda = hilda_VegFrac(site = site, lat = EddyCovariance_FLUXNET_data.tower_lat.values, lon= EddyCovariance_FLUXNET_data.tower_lon.values, config=exp_config).hilda_VegFrac_proc()
    
    #%% Merge all data
    ds = xr.merge([EddyCovariance_FLUXNET_data, SWC_TS_data, Climate_FLUXNET_data, Climate_ERA5_data, gapfill_climate_data, 
                   TEA_T_E, MCD43A_data, MxD11A_data, 
                   insitu_AGB_data, gloBbiomass_AGB_data, hansen_treecover_data, 
                   soilGrids_soilTexture_data, DIST_FOREST_FRAC, 
                   Frac_tree_insitu, Frac_grass_insitu, Frac_shrub_insitu, Frac_crop_insitu, Frac_savanna_insitu, Frac_wetland_insitu,
                   Frac_tree_hilda, Frac_grass_hilda, Frac_shrub_hilda, Frac_crop_hilda, Frac_other_hilda])
     
    #%% harmonize dataset and metadata to variables
    ds = data_harm(data = ds, config=exp_config).data_harm_proc()
                
    #%% Add global attributes
    ds = ds.assign_attrs(title = "FLUXNET data set for Sindbad",
                         FLUXNET_version = exp_config["FLUXNET_version"],
                         created_by='Simon Besnard',
                         contact = 'sbesnard@bgc-jena.mpg.de',
                         SITE_ID = site,
                         PFT = PFT_[0],
                         LAST_DISTURBANCE_DATE = str(LAST_DISTURBANCE_DATE[0]),
                         creation_date=datetime.now().strftime("%d-%m-%Y %H:%M"),                         
                         creation_script=CurrentScript, 
                         _FillValue = np.nan)
        
    #%% Export data
    if not os.path.exists(exp_config['OutPath']['nc_file'] + "/daily/"):
        os.makedirs(exp_config['OutPath']['nc_file'] + "/daily/")    
    ds.to_netcdf(exp_config['OutPath']['nc_file'] + "/daily/" + site + "." + exp_config['start_date'].split('-')[0] + "-" + str(int(exp_config['end_date'].split('-')[0]) -1) + '.nc', mode='w')
    if exp_config["temporal_resolution"]['monthly']:
        ds_monthtly = ds.resample(time='1M').mean()
        if not os.path.exists(exp_config['OutPath']['nc_file'] + "/monthly/"):
            os.makedirs(exp_config['OutPath']['nc_file'] + "/monthly/")    
        ds_monthtly.to_netcdf(exp_config['OutPath']['nc_file'] + "/monthly/" + site + "." + exp_config['start_date'].split('-')[0] + "-" + str(int(exp_config['end_date'].split('-')[0]) -1) + '.nc', mode='w')
    if exp_config["temporal_resolution"]['annual']:
       ds_annual =  ds.resample(time='1Y').mean()
       if not os.path.exists(exp_config['OutPath']['nc_file'] + "/annual/"):
           os.makedirs(exp_config['OutPath']['nc_file'] + "/annual/")        
       ds_annual.to_netcdf(exp_config['OutPath']['nc_file'] + "/annual/" + site + "." + exp_config['start_date'].split('-')[0] + "-" + str(int(exp_config['end_date'].split('-')[0]) -1) + '.nc', mode='w')
       
    #%% Create some diagnostic plots
    if not os.path.exists(exp_config['OutPath']['figs'] + site):
        os.makedirs(exp_config['OutPath']['figs'] + site) 
    for var_ in list(ds.keys()):
        dat_ = ds[var_]
        fig, ax = plt.subplots(1, 1, figsize=(16, 11), gridspec_kw={'wspace': 0.35, 'hspace': 0.5})
        ax.plot(ds.time.values.reshape(-1), dat_.values.reshape(-1))
        plt.savefig(exp_config['OutPath']['figs'] + site + "/" + var_ + "." + exp_config['start_date'].split('-')[0] + "-" + str(int(exp_config['end_date'].split('-')[0]) -1) + '.png', dpi=300)
    print(f'Successfull data processing for site {site}')    
    
    
