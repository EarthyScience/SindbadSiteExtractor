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
    def __init__(self, site_info, fill_src, fill_tar, qc_thres, data_dict,
                 config):
        self.site = site_info['site_ID']
        self.data_tar = data_dict[fill_tar].copy()
        self.data_src = data_dict[fill_src].copy()
        self.qc_thres = qc_thres
        self.info_src = config['dataset'][fill_src]
        self.info_tar = config['dataset'][fill_tar]
        self.fill_tar = fill_tar
        self.fill_src = fill_src
        self.temporal_resolution = config["temporal_resolution"]
        self.config = config

    def gapfill(self):
        data_gapfilled = []
        vars_in_src = list(
            set(self.data_src.variables) - set(self.data_src.coords))

        data_src, data_tar = xr.align(self.data_src,
                                      self.data_tar,
                                      join="outer")

        for var_src in vars_in_src:
            var_tar = var_src.replace(self.info_src['var_suffix'],
                                      self.info_tar['var_suffix'])
            if var_tar.endswith('_'):
                var_tar = var_tar[:-1]

            data_src_var = data_src[var_src]
            data_tar_var = data_tar[var_tar]
            qc_var = self.info_tar['variables'][var_tar]['QC']

            try:
                data_tar_gfw_tmp = data_tar_var.where(
                    data_tar[qc_var].fillna(0) > self.qc_thres)
                shut.log_and_print('GapFill', self.site, f'CLIFF',
                                   f'{var_tar}[{qc_var} < {self.qc_thres}]',
                                   var_src, self.temporal_resolution, self.config)
            except:
                data_tar_gfw_tmp = data_tar_var
                shut.log_and_print('GapFill', self.site, f'CLIFF',
                                   f'{var_tar}[{qc_var} not found!!]', var_src,
                                   self.temporal_resolution, self.config)

            data = xr.where(np.isfinite(data_tar_gfw_tmp), data_tar_gfw_tmp,
                            data_src_var)

            data = data.rename(var_tar + '_' + self.info_src['var_suffix'] +
                               '_gfld')

            data.attrs["units"] = 'units not set'
            shut.log_datavar_info(f'{self.fill_src} to {self.fill_tar}',
                                  self.site, self.temporal_resolution, var_src,
                                  var_tar, 'units_src = units_tar',
                                  'ensure source/target consistency',
                                  [np.nan, np.nan], data, None)

            data_gapfilled.append(data)

        data_gapfilled = xr.merge(data_gapfilled)
        return data_gapfilled


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print(
        'Provider at: ',
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Provider: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
