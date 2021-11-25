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

class soilgrids:
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
        soilGrids_soilTexture_df = pd.read_csv(self.cubepath)
        src_data = []
        for var_name in list(self.vars.keys()):
            src_name =  self.vars[var_name]['sourceVariableName']
            logger.info(f'------{self.site}: target: {var_name}, src: {src_name}: {self.temporal_resolution}-----')
            print (f'{self.site}: target: {var_name}, src: {src_name}: {self.temporal_resolution}')
            src_units = self.vars[var_name]['sourceVariableUnit']
            tar_units = self.vars[var_name]['variableUnit']

            soilTexture_var = soilGrids_soilTexture_df.filter(like = self.vars[var_name]['sourceVariableName'], axis=1).filter(regex='^((?!Sync).)*$').loc[soilGrids_soilTexture_df['SiteID'] == self.site]
            soilGrids_soilTexture_data = soilTexture_var.reindex(sorted(soilTexture_var.columns), axis=1).values.reshape(-1).astype('d')
            
            data = xr.DataArray(np.repeat(np.nan, len(soilTexture_var.columns)).reshape(len(soilTexture_var.columns), 1, 1), 
                                                        dims=['depth_soilGrids','latitude', 'longitude'],
                                                        coords={'depth_soilGrids': np.arange(len(soilTexture_var.columns))+1,
                                                                'latitude': [self.lat],
                                                                'longitude': [self.lon]}).to_dataset(name = var_name)
            if len(soilGrids_soilTexture_data) != 0:
                data[var_name].values = np.ones_like(data[var_name].values) * soilGrids_soilTexture_data.reshape(soilGrids_soilTexture_data.size, 1, 1)

            data.attrs["variable_name"]=var_name
            units_scalar = self.vars[var_name]['source2sindbadUnit']
            data = shut.set_units(data, src_name, src_units, tar_units, units_scalar)
            src_data.append(data)
            shut.log_site_info(self.dataset, self.site, self.temporal_resolution, src_name, var_name, self.vars[var_name]["variableUnit"], self.vars[var_name]["variableUnit"], self.vars[var_name]["sourceVariableUnit"], self.vars[var_name]["bounds"], data[var_name] , None)
        src_dataset =  xr.merge(src_data)
        return src_dataset
        

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
