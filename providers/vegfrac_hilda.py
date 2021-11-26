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

class vegfrac_hilda:
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
            logger.warning(f"vegfrac_hilda only provides daily or hourly (resampled) data but the resolution is {self.temporal_resolution}. The vegetation fraction data will not be included in {self.temporal_resolution} data for {self.site}.")
            return None

        date_ = np.arange(np.datetime64(self.start_date), np.datetime64(self.end_date) + np.timedelta64(1,'D') ,dtype=dtype_metadata)        

        src_data = []

        for var_name in list(self.vars.keys()):

            hilda_data_all = xr.open_dataset(self.vars[var_name]['data_path']).sel(time = slice(int(self.start_date[0:4]),int(self.end_date[0:4])))
            hilda_data = hilda_data_all.sel(latitude = self.lat, longitude = self.lon, method='nearest').LULC_states
            hilda_data = np.round(hilda_data, -1)
            
            src_name =  self.vars[var_name]['sourceVariableName']
            
            shut.log_and_print(self.site, self.vars[var_name]['sourceDataProductName'], var_name, src_name, self.temporal_resolution)

            data = xr.Dataset({var_name: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})

            if len(hilda_data) > 0:
                vegfrac = hilda_data.where(hilda_data == int(self.vars[var_name]['class']), 0)
                vegfrac[vegfrac > 0] = 1
                annual_vegfrac = []
                for year_ in hilda_data.time:
                    date_y = np.arange(str(int(year_.values)) + '-01' + '-01', str(int(year_.values) +1) + '-01' + '-01'  ,dtype=dtype_metadata)
                    annual_vegfrac.append(xr.Dataset({var_name: xr.DataArray(data = int(vegfrac.sel(time= int(year_.values)).values), dims = ['time'], coords = {'time': date_y})}))
                data = xr.merge(annual_vegfrac)

            data.attrs["variable_name"]=var_name

            data = shut.set_units(data, src_name, self.vars[var_name]['sourceVariableUnit'], self.vars[var_name]['variableUnit'], self.vars[var_name]['source2sindbadUnit'])

            src_data.append(data)
            
            shut.log_site_info(self.dataset, self.site, self.temporal_resolution, src_name, var_name, self.vars[var_name]["variableUnit"], self.vars[var_name]["variableUnit"], self.vars[var_name]["sourceVariableUnit"], self.vars[var_name]["bounds"], data[var_name] , None)

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