#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
from fluxcom.providers import eddy_covariance as ec
from fluxcom.transformers import hourly_to_daily
from fluxcom.core.variables import Variable
import xarray as xr
import utils.shared_utils as shut
import utils.shared_transformers as shtf

import logging
logger = logging.getLogger(__name__)
class bgi_FLUXCOM:
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
        src_prov    = ec.eddy_covariance.EddyProvider(cubepath = self.cubepath,version = self.version, site = self.site, NEE_partitioning_method = None)
        src_vars = [_var.name for _var in src_prov.variables]
        src_data = []
        for var_name in list(self.vars.keys()):
            src_name =  self.vars[var_name]['sourceVariableName']
            if src_name in src_vars:
                logger.info(f'------{self.site}: target: {var_name}, src: {src_name}: {self.temporal_resolution}-----')
                print (f'{self.site}: target: {var_name}, src: {src_name}: {self.temporal_resolution}')
                src_units = self.vars[var_name]['sourceVariableUnit']
                src_units_l = src_units.lower()
                tar_units = self.vars[var_name]['variableUnit']
                tar_units_l = tar_units.lower()
                src_partitioning = self.vars[var_name]['partitioning']
                transform = []
                if self.temporal_resolution == 'daily':
                    if self.vars[var_name]['isEnergy'] and src_units_l.startswith('w') and tar_units_l.startswith('mj'):
                        transform.append(shtf.eDaySum())
                    elif self.vars[var_name]['isCarbon'] and src_units_l.startswith('umol') and tar_units_l.startswith('gc'):
                        transform.append(shtf.cDaySum())
                    elif self.vars[var_name]['isWater']:
                        transform.append(shtf.DaySum())
                    elif 'DayTime' in var_name:
                        transform.append(shtf.DaytimeMean(src_prov.get_data(Variable('SW_IN_POT'))))
                    elif 'DayMin' in var_name:
                        transform.append(shtf.DayMin())
                    elif 'DayMax' in var_name:
                        transform.append(shtf.DayMax())
                    elif 'DaySum' in var_name:
                        transform.append(shtf.DaySum())
                    elif 'wDayMean' in var_name:
                        transform.append(shtf.wDayMean(src_prov.get_data(Variable(self.vars[var_name]['weightVar']))))
                    else:
                        transform = hourly_to_daily
                elif self.temporal_resolution == 'hourly':
                    if self.vars[var_name]['isEnergy'] and src_units_l.startswith('w') and tar_units_l.startswith('mj'):
                        transform.append(shtf.wm_2_to_mjh_1())
                    elif self.vars[var_name]['isCarbon'] and src_units_l.startswith('umol') and tar_units_l.startswith('gc'):
                        transform.append(shtf.umolm_2s_1_to_gCh_1())
                    else:
                        transform = []
                if len(transform) > 0:
                    src_prov.add_transform(transform)
                if src_name == 'P':
                    data = src_prov.get_data(Variable(src_name))
                else:
                    data = src_prov.get_data(Variable(src_name, units=src_units, partitioning=src_partitioning))
                if self.temporal_resolution == 'hourly':
                    data = shut.flatten_hour_to_time(data)
                data = data.rename(var_name)
                data.attrs["variable_name"]=var_name
                units_scalar = self.vars[var_name]['source2sindbadUnit']
                data = shut.set_units(data, src_name, src_units, tar_units, units_scalar)
                src_data.append(data)
                shut.log_site_info(self.dataset, self.site, self.temporal_resolution, src_name, var_name, data.attrs["units"], self.vars[var_name]["variableUnit"], self.vars[var_name]["sourceVariableUnit"], self.vars[var_name]["bounds"], data , src_prov.transforms, partitioning=src_partitioning)
                src_prov.transforms=[]
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