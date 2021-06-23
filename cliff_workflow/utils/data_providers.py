#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 10:21:27 2021

@author: simon
"""
from fluxcom.providers import eddy_covariance as ec
from fluxcom.variable import Variable
import xarray as xr
import pandas as pd
import numpy as np
import sys
sys.path.append('/Net/Groups/BGI/work_3/sindbad/sindbad-preproc/cliff_workflow/')
from utils.data_structure import data_structure_temporal_NoDepth

class EddyCovariance_FLUXNET:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.cubepath = config["DataPath"]["fluxcom_cube"]
        self.version = config["FLUXNET_version"]
        self.vars = config["variables"]["EddyCovariance_vars"]["FLUXNET_EddyCovariance"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]        
     
    def EddyCovariance_FLUXNET_proc(self):
        
        #%% Process data
        ec_prov    = ec.eddy_covariance.EddyProvider(cubepath   = self.cubepath,
                                                     NEE_partitioning_method = None,
                                                     version    = self.version, 
                                                     site       =self.site)
        ec_vars = []
        for k, v in self.vars.items():
            ec_vars.append(Variable(v['sourceVariableName'], partitioning = v['partitioning']))   
        ec_data = ec_prov.get_data(ec_vars)
        ec_data = ec_data.rename(**{v:v+'_FLUXNET' for v in list(set(ec_data.variables)-set(ec_data.coords))})
        
        #%% Apply unit correction and check value bounds
        for var_ in list(self.vars.keys()):
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            ec_data[var_] = ec_data[var_] * unit_scalar
            #if (np.nanmin(ec_data[var_]) < self.vars[var_]['bounds'][0]) | (np.nanmax(ec_data[var_]) > self.vars[var_]['bounds'][1]):
            #    print("Values are out of bounds:\n"+var_)
        
        return ec_data#.sel(time = slice(self.start_date, self.end_date))
    
class Climate_FLUXNET:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.cubepath = config["DataPath"]["fluxcom_cube"]
        self.version = config["FLUXNET_version"]
        self.vars = config["variables"]["climate_vars"]['FLUXNET_climate']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
        
    def Climate_FLUXNET_proc(self):
        
        #%% Process data
        Climate_FLUXNET_prov    = ec.eddy_covariance.EddyProvider(cubepath = self.cubepath,
                                                                  version        = self.version, 
                                                                  site           =self.site)
        Climate_FLUXNET_vars = []
        for var_ in list(self.vars.keys()) :
           Climate_FLUXNET_vars.append(Variable(self.vars[var_]['sourceVariableName']))   
        Climate_FLUXNET_data = Climate_FLUXNET_prov.get_data(Climate_FLUXNET_vars)
        
        #%% Rename data
        Climate_FLUXNET_data = Climate_FLUXNET_data.rename(**{v:v+'_FLUXNET' for v in list(set(Climate_FLUXNET_data.variables)-set(Climate_FLUXNET_data.coords))})
                
        #%% Apply unit correction and check value bounds
        for var_ in list(self.vars.keys()):
            unit_scalar = int(self.vars[var_]['source2sindbadUnit'])
            Climate_FLUXNET_data[var_] = Climate_FLUXNET_data[var_] * unit_scalar
            #if (np.nanmin(Climate_FLUXNET_data[var_]) < self.vars[var_]['bounds'][0]) | (np.nanmax(Climate_FLUXNET_data[var_]) > self.vars[var_]['bounds'][1]):
            #    print("Values are out of bounds:\n"+var_)
        
        return Climate_FLUXNET_data#.sel(time = slice(self.start_date, self.end_date))
    
class insitu_VegFrac:
    def __init__(self, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']        
        self.datapath = config["DataPath"]["PFT"]
        self.vars = config["variables"]["Veg_Frac_vars"]['Veg_Frac_insitu']
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]        
        
    def insitu_VegFrac_proc(self):
        
        #%% Process data     
        PFT_df = pd.read_csv(self.datapath)
        PFT_ = PFT_df.loc[PFT_df['Site_ID'] == self.site]
        if PFT_.size == 0:
            PFT_ = np.nan        
        return PFT_ 

class data_harm:
    def __init__(self, data, config):
        self.data = data
        self.config = config
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
        
    def data_harm_proc(self):
        
        #%Sort variable alphabeticatly
        var_ = sorted(list(self.data.keys()))
        self.data = self.data[var_]    
        
        #%% Add metadata
        for dat_type in self.config["variables"].keys():
                               
                for source_ in list(self.config["variables"][dat_type].keys()):
                    
                    for var_ in list(self.config["variables"][dat_type][source_].keys()):
                        self.data[var_] = self.data[var_].assign_attrs(bounds = str([np.nanmin(self.data[var_]), np.nanmax(self.data[var_])]),
                                                               units = self.config["variables"][dat_type][source_][var_]['variableUnit'],
                                                               short_name = self.config["variables"][dat_type][source_][var_]['nameShort'],
                                                               long_name = self.config["variables"][dat_type][source_][var_]['nameLong'],
                                                               sourceDataProductName = self.config["variables"][dat_type][source_][var_]['sourceDataProductName'],
                                                               publication = self.config["variables"][dat_type][source_][var_]['publication'],
                                                               dataPath = self.config["variables"][dat_type][source_][var_]['dataPath'])
                     
        return self.data.sel(time = slice(self.start_date, self.end_date))