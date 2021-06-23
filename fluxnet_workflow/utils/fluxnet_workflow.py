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
import sys
import numpy as np
sys.path.append('/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/fluxnet_workflow/')
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
            
        #%% Process downscaled climate data - ERA5
        Climate_ERA5_data = data_providers.Climate_ERA5(site_info = site_info, config=self.exp_settings).Climate_ERA5_proc()
                
        #%% Create gapfilled climate data
        gapfill_climate_data = data_providers.gapfill_climate(Climate_site_data = Climate_FLUXNET_data, Climate_downscaled_data = Climate_ERA5_data, config = self.exp_settings).gapfill_climate_proc()
        
        #%% Process soil water content and soil temperature data
        SWC_TS_data = data_providers.SWC_TS_FLUXNET(site_info = site_info, config=self.exp_settings).SWC_TS_FLUXNET_proc()    
            
        #%% Add T and E Nelson et al. data
        TEA_T_E = data_providers.EvapTrans_Nelson2018(site_info = site_info, config=self.exp_settings).EvapTrans_Nelson2018_proc()    
        
        #%% Process MODIS MCD43A data
        MCD43A_data = data_providers.MODIS_MCD43A(site_info = site_info, config=self.exp_settings).MODIS_MCD43A_proc()    
        
        #%% Process MODIS MxD11A data
        MxD11A_data = data_providers.MODIS_MxD11A(site_info = site_info, config=self.exp_settings).MODIS_MxD11A_proc()        
        
        #%% Process in-situ AGB data
        insitu_AGB_data = data_providers.insitu_AGB(site_info = site_info, config=self.exp_settings).insitu_AGB_proc()
        
        #%% Process globbiomass AGB data
        gloBbiomass_AGB_data = data_providers.gloBbiomass_AGB(site_info = site_info, config=self.exp_settings).gloBbiomass_AGB_proc()
        
        #%% Process Hansen tree cover data
        hansen_treecover_data = data_providers.hansen_treecover(site_info = site_info, config=self.exp_settings).hansen_treecover_proc() 
        
        #%% Process Soilgrids data
        soilGrids_soilTexture_data = data_providers.soilGrids_soilTexture(site_info = site_info, config=self.exp_settings).soilGrids_soilTexture_proc() 
            
        #%% Process age data
        DIST_FOREST_FRAC, LAST_DISTURBANCE_DATE = data_providers.forest_age(site_info = site_info, config=self.exp_settings).forest_age_proc()
        
        #%% Process insitu PFT variables
        PFT_, Frac_Veg_insitu = data_providers.insitu_VegFrac(site_info = site_info, config=self.exp_settings).insitu_VegFrac_proc()    
        
        #%% Process hilda+data
        Frac_Veg_hilda = data_providers.hilda_VegFrac(site_info = site_info, config=self.exp_settings).hilda_VegFrac_proc()
        
        #%% Merge all data
        ds = xr.merge([EddyCovariance_FLUXNET_data, SWC_TS_data, 
                    Climate_FLUXNET_data, Climate_ERA5_data, gapfill_climate_data, 
                    TEA_T_E, MCD43A_data, MxD11A_data, 
                    insitu_AGB_data, gloBbiomass_AGB_data, hansen_treecover_data, 
                    soilGrids_soilTexture_data, DIST_FOREST_FRAC, 
                    Frac_Veg_insitu, Frac_Veg_hilda])
        
        #%% harmonize dataset and metadata to variables
        ds = data_providers.data_harm(data = ds, config=self.exp_settings).data_harm_proc()
        
        #%% Aggregate to monthly and annual scale
        ds_monthtly = ds.resample(time='1M').mean()
        ds_annual =  ds.resample(time='1Y').mean()            
        
        #%% Transform time dim to days since time dimension 
        date_ = np.arange(np.datetime64('1592-10-15'), np.datetime64('2019-12-31') + np.timedelta64(1,'D') ,dtype="M8[D]")
        days_since = np.arange(0, len(date_))
        days_since = days_since[np.in1d(date_, ds.time)]
        ds['time'] = days_since
        ds['time'] = ds['time'].assign_attrs(units = 'days since 1592-10-15 00:00:00', 
                                             calendar = 'proleptic_gregorian') 
                    
        #%% Add global attributes
        CurrentScript = os.path.basename('/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/fluxnet_workflow/fluxnet_workflow.py')
        ds = ds.assign_attrs(title = "FLUXNET data set for Sindbad",
                            FLUXNET_version = self.exp_settings["FLUXNET_version"],
                            created_by='Simon Besnard',
                            contact = 'sbesnard@bgc-jena.mpg.de',
                            SITE_ID = site_info["site_ID"],
                            latitude = site_info["latitude"],                         
                            longitude = site_info["longitude"],
                            start_year =  self.exp_settings['start_date'].split('-')[0],
                            end_year =  self.exp_settings['end_date'].split('-')[0],
                            PFT = PFT_[0],
                            LAST_DISTURBANCE_DATE = str(LAST_DISTURBANCE_DATE[0]),                            
                            creation_date=datetime.now().strftime("%d-%m-%Y %H:%M"),                         
                            creation_script=CurrentScript, 
                            _FillValue = np.nan)
        
        #%% Export Metadata
        #dict_data = ds.to_dict(data=False)
        #dict_data = {site: {'start_year': self.exp_settings['start_date'].split('-')[0], 
        #        'end_year': self.exp_settings['end_date'].split('-')[0],
        #        'latitude': dict_data['attrs']['latitude'],
        #        'longitude': dict_data['attrs']['longitude'],
        #        'PFT': dict_data['attrs']['PFT']}}     
        #with open(self.exp_settings['OutPath']['nc_file'] + "/daily/" + 'site_info.json', 'a') as fp:
        #    json.dump(dict_data, fp, indent=2)
        #    fp.write(",\n")
        
        #%% Export data
        if not os.path.exists(self.exp_settings['OutPath']['nc_file'] + "/daily/"):
            os.makedirs(self.exp_settings['OutPath']['nc_file'] + "/daily/")    
        ds.to_netcdf(self.exp_settings['OutPath']['nc_file'] + "/daily/" + site + "." + self.exp_settings['start_date'].split('-')[0] + "-" + str(int(self.exp_settings['end_date'].split('-')[0])) + '.nc', mode='w')
        if not os.path.exists(self.exp_settings['OutPath']['nc_file'] + "/monthly/"):
            os.makedirs(self.exp_settings['OutPath']['nc_file'] + "/monthly/")    
        ds_monthtly.to_netcdf(self.exp_settings['OutPath']['nc_file'] + "/monthly/" + site + "." + self.exp_settings['start_date'].split('-')[0] + "-" + str(int(self.exp_settings['end_date'].split('-')[0])) + '.nc', mode='w')
        if not os.path.exists(self.exp_settings['OutPath']['nc_file'] + "/annual/"):
            os.makedirs(self.exp_settings['OutPath']['nc_file'] + "/annual/")        
        ds_annual.to_netcdf(self.exp_settings['OutPath']['nc_file'] + "/annual/" + site + "." + self.exp_settings['start_date'].split('-')[0] + "-" + str(int(self.exp_settings['end_date'].split('-')[0])) + '.nc', mode='w')
        
        #%% Create diagnostic plots
        if not os.path.exists(self.exp_settings['OutPath']['figs'] + site):
            os.makedirs(self.exp_settings['OutPath']['figs'] + site) 
        for var_ in list(ds.keys()):
            dat_ = ds[var_]
            if ('depth_FLUXNET' in dat_.dims) & ('time' in dat_.dims):
                for layer_ in dat_.depth_FLUXNET:        
                    fig, ax = plt.subplots(1, 1, figsize=(16, 11), gridspec_kw={'wspace': 0.35, 'hspace': 0.5})
                    ax.plot(dat_.sel(depth_FLUXNET = layer_).time.values.reshape(-1), dat_.sel(depth_FLUXNET = layer_).values.reshape(-1))
                    ax.set_xlabel('time', size=18)
                    ax.set_ylabel(var_ + '_' + str(layer_.values) + ' [' + dat_.units + ']', size=18)        
                    ax.tick_params(labelsize=16)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    plt.savefig(self.exp_settings['OutPath']['figs'] + site + "/" + var_ + '_' + str(layer_.values) + "." + self.exp_settings['start_date'].split('-')[0] + "-" + str(int(self.exp_settings['end_date'].split('-')[0])) + '.png', dpi=300)
            elif ('time' in dat_.dims) & ('depth_FLUXNET' not in dat_.dims):
                fig, ax = plt.subplots(1, 1, figsize=(16, 11), gridspec_kw={'wspace': 0.35, 'hspace': 0.5})
                ax.plot(dat_.time.values.reshape(-1), dat_.values.reshape(-1))
                ax.set_xlabel('time', size=18)
                ax.set_ylabel(var_ + ' [' + dat_.units + ']', size=18)        
                ax.tick_params(labelsize=16)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.savefig(self.exp_settings['OutPath']['figs'] + site + "/" + var_ + "." + self.exp_settings['start_date'].split('-')[0] + "-" + str(int(self.exp_settings['end_date'].split('-')[0])) + '.png', dpi=300)
            elif ('time' not in dat_.dims) & ('depth_soilGrids' in dat_.dims):
                for layer_ in dat_.depth_soilGrids:
                    fig, ax = plt.subplots(1, 1, figsize=(16, 11), gridspec_kw={'wspace': 0.35, 'hspace': 0.5})
                    ax.scatter(var_ + '_' + str(layer_.values), 
                        dat_.sel(depth_soilGrids = layer_).values.reshape(-1), s =1000)
                    ax.set_ylabel(var_ + '_' + str(layer_.values) + ' [' + dat_.units + ']', size=18)        
                    ax.tick_params(labelsize=16)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    plt.savefig(self.exp_settings['OutPath']['figs'] + site + "/" + var_ + '_' + str(layer_.values) + "." + self.exp_settings['start_date'].split('-')[0] + "-" + str(int(self.exp_settings['end_date'].split('-')[0])) + '.png', dpi=300)
        print(f'Successfull data processing for site {site}')    
    
    
