#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: @dr-ko
"""
from fluxcom.core.variables import Variable
from fluxcom.providers import MCD43A


from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)

import numpy as np
import xarray as xr

def extract(dataset, site_info, config):

    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='daily'):
        return None

    src_prov    = MCD43A(cubepath=bxtr.flx_cubepath, site=bxtr.site)
    src_vars = [_var.name for _var in src_prov.variables]
    src_data = []
    src_data_dict = {}
    for tar_name in bxtr.vars_list:
        src_name =  bxtr.vars[tar_name]['sourceVariableName']
        if src_name in src_vars:

            data = src_prov.get_data(Variable(src_name, units=bxtr.vars[tar_name]['sourceVariableUnit'], partitioning=bxtr.vars[tar_name]['partitioning']))

            data = data.rename(tar_name)
            
            data = bxtr.convert_units(data, tar_name)

            src_data_dict[tar_name] = data

    # merge the data manually (rough method but works.. be careful of hardcoded variable names. things may break if there are changes in extractor's json settings)

    data = xr.where(np.isfinite(src_data_dict['mergedNDVI_MCD43A']), src_data_dict['kNDVI_MCD43A'], src_data_dict['NDVI_MCD43A'])
    tar_name = 'mergedNDVI_MCD43A'
    bxtr.log_var_start(tar_name)
    # needs renaming and resetting the attributes because xarray likes deleting things
    data = data.rename(tar_name)
    data = bxtr.convert_units(data, tar_name)
    bxtr.log_var_end(data, tar_name, src_prov.transforms)
    src_data_dict[tar_name] = data.copy()

    src_dataset =  bxtr.merge_and_format([src_data_dict[tar_name]])   
    return src_dataset


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)