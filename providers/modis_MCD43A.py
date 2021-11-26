#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
from fluxcom.core.variables import Variable
from fluxcom.providers.modis.refl.modis_MCD43A import MCD43A
import utils.shared_utils as shut
import xarray as xr

import logging
logger = logging.getLogger(__name__)

class modis_MCD43A:
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
        if self.temporal_resolution != 'daily':
            logger.warning(f"modis_MCD43A only provides daily data but the resolution is {self.temporal_resolution}. The MCD43A reflectances will not be included in {self.temporal_resolution} data for {self.site}.")
            return None
        src_prov    = MCD43A(cubepath=self.flx_cubepath, site=self.site)
        src_vars = [_var.name for _var in src_prov.variables]
        src_data = []
        for var_name in list(self.vars.keys()):
            src_name =  self.vars[var_name]['sourceVariableName']
            if src_name in src_vars:

                src_partitioning = self.vars[var_name]['partitioning']
                shut.log_and_print(self.site, self.vars[var_name]['sourceDataProductName'], var_name, src_name, self.temporal_resolution)

                data = src_prov.get_data(Variable(src_name, units=self.vars[var_name]['sourceVariableUnit'], partitioning=src_partitioning))
                data = data.rename(var_name)
                data.attrs["variable_name"]=var_name
                
                data = shut.set_units(data, src_name, self.vars[var_name]['sourceVariableUnit'], self.vars[var_name]['variableUnit'], self.vars[var_name]['source2sindbadUnit'])
                src_data.append(data)

                shut.log_site_info(self.dataset, self.site, self.temporal_resolution, src_name, var_name, data.attrs["units"], self.vars[var_name]["variableUnit"], self.vars[var_name]["sourceVariableUnit"], self.vars[var_name]["bounds"], data , src_prov.transforms, partitioning=src_partitioning)


        src_dataset =   xr.merge(src_data)
        src_dataset =   shut.data_structure_temporal_NoDepth(src_dataset, self.lat, self.lon)                  
        return src_dataset


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)