#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
from numpy.core.arrayprint import format_float_positional
import xarray as xr
import numpy as np
import utils.shared_utils as shut

import logging
logger = logging.getLogger(__name__)

class cliff_gapfill:
    def __init__(self, site_info, fill_src, fill_tar, qc_thres, data_dict,config):
        self.site = site_info['site_ID']
        self.data_tar = data_dict[fill_tar].copy()
        self.data_src = data_dict[fill_src].copy()
        self.qc_thres = qc_thres
        self.info_src = config['dataset'][fill_src]  
        self.info_tar = config['dataset'][fill_tar]  
        self.fill_tar = fill_tar
        self.fill_src = fill_src
        self.temporal_resolution = config["temporal_resolution"]

    def process(self):
        data_gapfilled = []
        vars_in_src = list(set(self.data_src.variables)-set(self.data_src.coords))
        data_src, data_tar = xr.align(self.data_src, self.data_tar, join="outer")
        logger.info(f'{self.site}: cliff_gapfill: target: {self.fill_tar}, source: {self.fill_src} using quality flag threshold of {self.qc_thres}')
        for var_src in vars_in_src:
            var_tar = var_src.replace(self.info_src['var_suffix'],self.info_tar['var_suffix'])
            if var_tar.endswith('_'):
                var_tar = var_tar[:-1]
            data_src_var = data_src[var_src]
            data_tar_var = data_tar[var_tar]
            try:
                qc_var = var_tar+'_QC'
                if len(self.info_tar['var_suffix']) > 0:
                    qc_var = qc_var + '_' + self.info_tar['var_suffix']
                data_tar_gfw_tmp = data_tar_var.where(data_tar[qc_var] > self.qc_thres)
            except:
                data_tar_gfw_tmp = data_tar_var
                logger.info(f"No QC variable {qc_var} found for {var_tar}. QC masking will not be applied")
            
            data = xr.where(np.isfinite(data_tar_gfw_tmp), data_tar_gfw_tmp, data_src_var)
            data = data.rename(var_tar+'_'+self.info_src['var_suffix']+'_gfld')
            shut.log_site_info(f'{self.fill_src} to {self.fill_tar}', self.site, self.temporal_resolution, var_src, var_tar, 'units assumed same', 'units assumed same', 'units assumed same', [np.nan,np.nan], data , None)

            data_gapfilled.append(data)

        data_gapfilled =   xr.merge(data_gapfilled)
        return data_gapfilled

        

if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print('Provider at: ',os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))) 
    print ('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
