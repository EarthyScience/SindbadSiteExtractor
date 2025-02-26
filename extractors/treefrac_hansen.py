#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: dr-ko
"""
import pandas as pd
import xarray as xr
import numpy as np

from extractors.BasexTractor import BasexTractor
import logging
logger = logging.getLogger(__name__)


def extract(dataset, site_info, config):
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    date_ = bxtr.get_date_vec()

    src_data = []
    for tar_name in bxtr.vars_list:

        src_df = pd.read_csv(bxtr.vars[tar_name]['data_path'])
        src_name =  bxtr.vars[tar_name]['sourceVariableName']
        sel_data = np.array(src_df.loc[src_df['SiteID'] == bxtr.site][src_name], dtype=float)
        
        bxtr.log_var_start(tar_name)
        
        data = xr.Dataset({tar_name: xr.DataArray(data = np.nan, dims = ['time'], coords = {'time': date_})})

        if len(sel_data) == 0:
            data[tar_name].values = np.zeros_like(data[tar_name].values)
        else:
            data[tar_name].values = np.ones_like(data[tar_name].values) * sel_data
        
        data = bxtr.convert_units(data, tar_name)
        bxtr.log_var_end(data, tar_name, None)
        src_data.append(data)

    src_dataset =  bxtr.merge_and_format(src_data)   
                  
    return src_dataset
        

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
