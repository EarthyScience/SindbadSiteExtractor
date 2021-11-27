#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
from fluxcom.providers import eddy_covariance as ec
from fluxcom.core.variables import Variable
import xarray as xr
import numpy as np
import utils.shared_utils as shut

from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)


def extract(dataset, site_info, config):

    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    src_prov    = ec.eddy_covariance.EddyProvider(cubepath = bxtr.flx_cubepath,version = bxtr.version, site = bxtr.site, NEE_partitioning_method = None)
    src_vars = [_var.name for _var in src_prov.variables]
    src_data = []
    variants = 'shallow deep'.split()
    for tar_name in list(bxtr.vars.keys()):
        src_name =  bxtr.vars[tar_name]['sourceVariableName']
        if src_name in src_vars:
            
            bxtr.log_var_start(tar_name)
            
            transform = bxtr.get_transform(tar_name,src_prov=src_prov)
            src_prov.add_transform(transform)
            shut.log_and_print('xTrct', bxtr.site, bxtr.vars[tar_name]['sourceDataProductName'], tar_name, src_name, bxtr.temporal_resolution)
            
            
            data = []
            for depth in variants:
                data_tmp = src_prov.get_data(Variable(src_name, units=bxtr.vars[tar_name]['sourceVariableUnit'], partitioning=bxtr.vars[tar_name]['partitioning'], depth=depth))
                data.append(data_tmp)

            if bxtr.temporal_resolution == 'hourly':
                data = [shut.flatten_hour_to_time(_data) for _data in data]

            data = xr.merge(data)
            data = xr.DataArray(data.to_array().values.reshape(len(list(data.keys())),-1,1,1), 
                        dims=['depth_FLUXNET', 'time', 'latitude', 'longitude'],
                        coords={'depth_FLUXNET': np.arange(len(list(data.keys()))) +1,
                                'time': data.time.values,
                                'latitude': [bxtr.lat],
                                'longitude': [bxtr.lon]})
                                 
            data = data.rename(tar_name)

            data = bxtr.convert_units(data, tar_name)
            src_data.append(data)

            bxtr.log_var_end(data, tar_name, transform)
            src_prov.transforms = []
    src_dataset =  bxtr.merge_and_format(src_data)   
    return src_dataset
                    

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)