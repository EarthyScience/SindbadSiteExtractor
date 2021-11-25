#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
from fluxcom.core.variables import Variable
from fluxcom.providers.modis.LST.modis_MxD11A1 import MxD11A1
import utils.shared_utils as shut
import utils.shared_transformers as shtf
import xarray as xr
import numpy as np
import logging
logger = logging.getLogger(__name__)

class modis_MxD11A1:
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
        if self.temporal_resolution != 'daily':
            logger.warning(f"modis_MxD11A1 only provides daily data but the resolution is {self.temporal_resolution}. The MxD11A1 LST will not be included in {self.temporal_resolution} data for {self.site}.")
            return []
        src_prov    = MxD11A1(cubepath=self.cubepath, site=self.site, transforms=[shtf.k_to_degC()])

        src_vars = [_var.name for _var in src_prov.variables]
        src_data = []
        variants = 'Day Night'.split()
        for var_name in list(self.vars.keys()):
            src_name =  self.vars[var_name]['sourceVariableName']
            if src_name in src_vars:
                src_units = self.vars[var_name]['sourceVariableUnit']
                tar_units = self.vars[var_name]['variableUnit']
                satellite = self.vars[var_name]['satellite']
                data_tmp = {}
                for daynight in variants:
                    print (f'{self.site}: target: {var_name}, src: {src_name}, satellite: {satellite}, daynight: {daynight}:: {self.temporal_resolution}')
                    data_tmp[daynight] = src_prov.get_data(Variable(src_name, units=src_units, day_night=daynight, satellite=satellite))  - 273.15
                if 'DayTime' in var_name:
                    data = data_tmp['Day']
                else:
                    data = xr.merge([data_tmp['Day'], data_tmp['Night']]).to_array(dim='new').mean('new')


                src_partitioning = self.vars[var_name]['partitioning']

                data = data.rename(f'{var_name}')
                data.attrs["variable_name"]=var_name
                data.attrs["units"] = "deg C"
                
                units_scalar = self.vars[var_name]['source2sindbadUnit']
                data = shut.set_units(data, src_name, src_units, tar_units, units_scalar)
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