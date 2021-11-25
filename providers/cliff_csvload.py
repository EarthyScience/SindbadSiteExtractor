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
import utils.setup_config as set_conf
import os
import logging
logger = logging.getLogger(__name__)

class cliff_csvload:
    def __init__(self, dataset, site_info, config):
        self.site = site_info['site_ID']
        self.lat = site_info['latitude']
        self.lon = site_info['longitude']
        self.dataset = dataset
        self.cubepath = config['dataset'][dataset]["cube_data_path"]
        self.version = config["FLUXNET_version"]
        self.vars = config["dataset"][dataset]["variables"]
        self.source = config["dataset"][dataset]["src_cliff"]
        self.src_start_date = config["dataset"][dataset]["src_start_date"]
        self.src_end_date = config["dataset"][dataset]["src_end_date"]
        self.start_date = config["start_date"]
        self.end_date = config["end_date"]
        self.temporal_resolution = config["temporal_resolution"]
        
    def process(self):
        in_sub_path = set_conf.get_cliff_dirname(self.version, self.temporal_resolution)
        site_data_path = os.path.join(self.cubepath, in_sub_path, self.source, f'{self.site}.{self.src_start_date.split("-")[0]}.{self.src_end_date.split("-")[0]}.{self.temporal_resolution}.csv')
        if self.temporal_resolution == 'daily':
            dtype_metadata = "M8[D]"
        elif self.temporal_resolution == 'hourly':
            dtype_metadata = "M8[h]"
        else:
            logger.warning(f"cliff_csvload only provides daily or hourly (resampled) data but the resolution is {self.temporal_resolution}. The CLIFF variables will not be included in {self.temporal_resolution} data for {self.site}.")
            return []

        src_df = pd.read_csv(site_data_path)
            
        date_ = np.arange(np.datetime64(self.src_start_date), np.datetime64(self.src_end_date) + np.timedelta64(1,'D') ,dtype=dtype_metadata)        

        src_data = []

        logger.info (f'------{self.site}: cliff_csvload:: getting downscaled data  using \n\t {site_data_path}-----')
        for var_name in list(self.vars.keys()):
            src_name =  self.vars[var_name]['sourceVariableName']
            for ext in ['_DayTime', '_DayMean', '_DayMin' , '_DayMax']:
                if ext in var_name:
                    src_name = src_name + ext
            src_units = self.vars[var_name]['sourceVariableUnit']
            tar_units = self.vars[var_name]['variableUnit']
            if src_name in src_df:
                sel_data = src_df[src_name].values
                data = xr.Dataset({var_name: xr.DataArray(data = sel_data, dims = ['time'], coords = {'time': date_})})

                
                # data.attrs["variable_name"]=var_name
                units_scalar = self.vars[var_name]['source2sindbadUnit']
                data = shut.set_units(data[var_name], src_name, src_units, tar_units, units_scalar)
                src_data.append(data)

                shut.log_site_info(self.dataset, self.site, self.temporal_resolution, src_name, var_name, self.vars[var_name]["variableUnit"], self.vars[var_name]["variableUnit"], self.vars[var_name]["sourceVariableUnit"], self.vars[var_name]["bounds"], data , None)
            else:
                logger.warning(f'{self.site}: target: {var_name}, src: {src_name} was not found in the csv data frame in {site_data_path}')

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