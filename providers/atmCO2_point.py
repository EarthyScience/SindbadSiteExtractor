#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
import xarray as xr
import numpy as np
import utils.shared_utils as shut

import logging
logger = logging.getLogger(__name__)

class atmCO2_point:
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
        src_data = []
        for var_name in list(self.vars.keys()):
            src_name =  self.vars[var_name]['sourceVariableName']
            logger.info(f'------{self.site}: target: {var_name}, src: {src_name}: {self.temporal_resolution}-----')
            print (f'{self.site}: target: {var_name}, src: {src_name}: {self.temporal_resolution}')
            src_units = self.vars[var_name]['sourceVariableUnit']
            tar_units = self.vars[var_name]['variableUnit']

            data = xr.open_dataset(self.vars[var_name]['data_path'])
            #%% Gapfill data
            data =  data.where(np.isfinite(data)).interpolate_na(dim='time',method='nearest').compute()
            data =  data.rename(**{src_name: var_name})
            units_scalar = self.vars[var_name]['source2sindbadUnit']
            data = shut.set_units(data[var_name], src_name, src_units, tar_units, units_scalar)
            src_data.append(data)
            shut.log_site_info(self.dataset, self.site, self.temporal_resolution, src_name, var_name, self.vars[var_name]["variableUnit"], self.vars[var_name]["variableUnit"], self.vars[var_name]["sourceVariableUnit"], self.vars[var_name]["bounds"], data , None)
        src_dataset =  xr.merge(src_data)
        src_dataset = shut.data_structure_temporal_NoDepth(src_dataset, self.lat, self.lon)                
        return src_dataset
        

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
