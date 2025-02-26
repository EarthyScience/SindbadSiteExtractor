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
    
    src_data = []
    for tar_name in bxtr.vars_list:
        
        bxtr.log_var_start(tar_name)

        src_df = pd.read_csv(bxtr.vars[tar_name]['data_path'])
        soilTexture_var = src_df.filter(like = bxtr.vars[tar_name]['sourceVariableName'], axis=1).filter(regex='^((?!Sync).)*$').loc[src_df['SiteID'] == bxtr.site]
        soilGrids_soilTexture_data = soilTexture_var.reindex(sorted(soilTexture_var.columns), axis=1).values.reshape(-1).astype('d')

        
        data = xr.DataArray(np.ones((len(soilTexture_var.columns), 1, 1)) *np.nan, dims=['depth_soilGrids','latitude', 'longitude'], coords={'depth_soilGrids': np.arange(len(soilTexture_var.columns))+1,'latitude': [bxtr.lat],'longitude': [bxtr.lon]}).to_dataset(name = tar_name)

        if len(soilGrids_soilTexture_data) != 0:
            data[tar_name].values = np.ones_like(data[tar_name].values) * soilGrids_soilTexture_data.reshape(soilGrids_soilTexture_data.size, 1, 1)
        else:
            logger.warning(f"::MISSING:: variable {tar_name} has no data in source {bxtr.vars[tar_name]['data_path']} for {bxtr.site}. NaN values will be set.")


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
