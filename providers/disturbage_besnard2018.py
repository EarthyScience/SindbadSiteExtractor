#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A data provicer for forest age and disturbance
Created on Nov 23 2021

@author: sujan
"""
import sys, os
sys.path.append(os.path.join(os.getcwd(), '../'))
import pandas as pd
import xarray as xr
import numpy as np
import utils.shared_utils as shut

import logging
logger = logging.getLogger(__name__)

class disturbage_besnard2018:
    def __init__(self, dataset, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']
        self.dataset = dataset
        self.flx_cubepath = config["fluxcom_cube_path"]
        self.version = config["FLUXNET_version"]
        self.vars = config["dataset"][dataset]["variables"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
        self.temporal_resolution = config["temporal_resolution"]
        
    def process(self):

        if self.temporal_resolution == 'daily':
            dtype_metadata = "M8[D]"
        elif self.temporal_resolution == 'hourly':
            dtype_metadata = "M8[h]"
        else:
            logger.warning(f"disturbage_besnard2018 only provides daily or hourly (resampled) data but the resolution is {self.temporal_resolution}. The vegetation fraction data will not be included in {self.temporal_resolution} data for {self.site}.")
            return None
        date_ = np.arange(np.datetime64(self.start_date), np.datetime64(self.end_date) + np.timedelta64(1,'D') ,dtype=dtype_metadata)        


        src_data = []

        for var_name in list(self.vars.keys()):
            src_df = pd.read_csv(self.vars[var_name]['data_path'])
            last_disturbance_on =  src_df.loc[src_df['Site_ID'] == self.site]['Plantation_Date_max'].values.astype(np.datetime64)
            
            src_name =  self.vars[var_name]['sourceVariableName']

            shut.log_and_print(self.site, self.vars[var_name]['sourceDataProductName'], var_name, src_name, self.temporal_resolution)

            data = xr.Dataset({var_name: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})

            if len(last_disturbance_on) == 0:
                last_disturbance_on = ['undisturbed']
            else:
                forest_age_data = np.concatenate([(day_ - last_disturbance_on).astype('int') / 365.25 for day_ in date_] )
                forest_age_data = xr.DataArray(forest_age_data, dims=['time'],
                                        coords={'time': date_})
                forest_age_data = forest_age_data.where(forest_age_data==0)
                forest_age_data = forest_age_data.where(np.isfinite(forest_age_data),1)
                forest_age_data.values = 1-forest_age_data.values              

                data[var_name] = forest_age_data
                
            data.attrs["variable_name"]=var_name

            data = shut.set_units(data, src_name, self.vars[var_name]['sourceVariableUnit'], self.vars[var_name]['variableUnit'], self.vars[var_name]['source2sindbadUnit'])

            src_data.append(data)
            
            shut.log_site_info(self.dataset, self.site, self.temporal_resolution, src_name, var_name, self.vars[var_name]["variableUnit"], self.vars[var_name]["variableUnit"], self.vars[var_name]["sourceVariableUnit"], self.vars[var_name]["bounds"], data[var_name] , None)

        src_dataset =  xr.merge(src_data)
        src_dataset =   shut.data_structure_temporal_NoDepth(src_dataset, self.lat, self.lon)                  
        src_dataset = src_dataset.assign_attrs(last_disturbance_on= last_disturbance_on[0])
        return src_dataset
        

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)  