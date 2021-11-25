#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
import pandas as pd
import xarray as xr
import numpy as np
import utils.shared_utils as shut

import logging
logger = logging.getLogger(__name__)

class treefrac_hansen:
    def __init__(self, dataset, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']
        self.dataset = dataset
        self.cubepath = config['dataset'][dataset]["cube_data_path"]
        self.version = config["FLUXNET_version"]
        self.vars = config["dataset"][dataset]["variables"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
        self.temporal_resolution = config["temporal_resolution"]
        
    def process(self):
        src_df = pd.read_csv(self.cubepath)
        
        if self.temporal_resolution == 'daily':
            dtype_metadata = "M8[D]"
        elif self.temporal_resolution == 'hourly':
            dtype_metadata = "M8[h]"
        else:
            logger.warning(f"treefrac_hansen only provides daily or hourly (resampled) data but the resolution is {self.temporal_resolution}. The vegetation fraction data will not be included in {self.temporal_resolution} data for {self.site}.")
            return []
            
        date_ = np.arange(np.datetime64(self.start_date), np.datetime64(self.end_date) + np.timedelta64(1,'D') ,dtype=dtype_metadata)        

        src_data = []

        for var_name in list(self.vars.keys()):
            src_name =  self.vars[var_name]['sourceVariableName']
            sel_data = np.array(src_df.loc[src_df['SiteID'] == self.site][src_name], dtype=float)
            logger.info(f'------{self.site}: target: {var_name}, src: {src_name}: {self.temporal_resolution}-----')
            print (f'{self.site}: target: {var_name}, src: {src_name}: {self.temporal_resolution}')
            src_units = self.vars[var_name]['sourceVariableUnit']
            tar_units = self.vars[var_name]['variableUnit']
            data = xr.Dataset({var_name: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})
            if len(sel_data) == 0:
                data[var_name].values = np.zeros_like(data[var_name].values)
            else:
                data[var_name].values = np.ones_like(data[var_name].values) * sel_data
            
            # data.attrs["variable_name"]=var_name
            units_scalar = self.vars[var_name]['source2sindbadUnit']
            data = shut.set_units(data[var_name], src_name, src_units, tar_units, units_scalar)
            src_data.append(data)

            shut.log_site_info(self.dataset, self.site, self.temporal_resolution, src_name, var_name, self.vars[var_name]["variableUnit"], self.vars[var_name]["variableUnit"], self.vars[var_name]["sourceVariableUnit"], self.vars[var_name]["bounds"], data , None)
        src_dataset =  xr.merge(src_data)
        src_dataset =   shut.data_structure_temporal_NoDepth(src_dataset, self.lat, self.lon)                  
        return src_dataset
        

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
