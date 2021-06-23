#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 14:49:25 2019

@author: simon
"""
#%% Load library
import xarray as xr
import os
from datetime import datetime
import sys
import numpy as np
sys.path.append('/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/cliff_workflow/')
from utils import data_providers

class DataProcessor:    
    def __init__(self, exp_settings=None):
        self.exp_settings = exp_settings
        
    #%% Run data processing for each site
    def compile_data(self, site):
        
        print(f'Processing data for site {site}')
        #%% Get site information
        lat = xr.open_zarr(self.exp_settings["DataPath"]['fluxcom_cube']).site_var__v000__tower_lat.sel(site = site).values
        lon = xr.open_zarr(self.exp_settings["DataPath"]['fluxcom_cube']).site_var__v000__tower_lon.sel(site = site).values    
        site_info = {"site_ID": site, "latitude": lat, "longitude": lon}    
        
        #%% Process eddy-covariance data
        EddyCovariance_FLUXNET_data = data_providers.EddyCovariance_FLUXNET(site_info= site_info, config=self.exp_settings).EddyCovariance_FLUXNET_proc()
        
        #%% Process FLUXNET climate data
        Climate_FLUXNET_data = data_providers.Climate_FLUXNET(site_info= site_info, config=self.exp_settings).Climate_FLUXNET_proc()
        
        #%% Merge all data
        ds = xr.merge([EddyCovariance_FLUXNET_data,Climate_FLUXNET_data])
        
        #%% harmonize dataset and metadata to variables
        ds = data_providers.data_harm(data = ds, config=self.exp_settings).data_harm_proc()
                            
        #%% Add global attributes
        CurrentScript = os.path.basename('/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/fluxnet_workflow/fluxnet_workflow.py')
        ds = ds.assign_attrs(title = "Hourly FLUXNET data set for Cliff downscaling",
                            FLUXNET_version = self.exp_settings["FLUXNET_version"],
                            created_by='Simon Besnard',
                            contact = 'sbesnard@bgc-jena.mpg.de',
                            SITE_ID = site_info["site_ID"],
                            latitude = site_info["latitude"],                         
                            longitude = site_info["longitude"],
                            start_year =  self.exp_settings['start_date'].split('-')[0],
                            end_year =  self.exp_settings['end_date'].split('-')[0],
                            creation_date=datetime.now().strftime("%d-%m-%Y %H:%M"),                         
                            creation_script=CurrentScript, 
                            _FillValue = np.nan)
        
        #%% Export data
        if not os.path.exists(self.exp_settings['OutPath']['nc_file'] + "/hourly/" + self.exp_settings['FLUXNET_version']):
            os.makedirs(self.exp_settings['OutPath']['nc_file'] + "/hourly/" + self.exp_settings['FLUXNET_version'])    
        ds.to_netcdf(self.exp_settings['OutPath']['nc_file'] + "/hourly/" + self.exp_settings['FLUXNET_version'] + '/' + site + "." + self.exp_settings['start_date'].split('-')[0] + "-" + str(int(self.exp_settings['end_date'].split('-')[0])) + '.nc', mode='w')
        print(f'Successfull data processing for site {site}')    
    
    
