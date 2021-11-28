#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Nov 23 2021

@author: sujan
"""
import pandas as pd
import xarray as xr
import utils.setup_config as set_conf
import os
import logging
logger = logging.getLogger(__name__)
from extractors.BasexTractor import BasexTractor
import sys

def extract(dataset, site_info, config):
    bxtr = BasexTractor(dataset=dataset, site_info=site_info, config=config)

    if not bxtr.is_resolution_supported(provided_reso='hourly'):
        return None

    date_ = bxtr.get_date_vec()
    src_data = []

    for tar_name in bxtr.vars_list:
        site_data_path = os.path.join(
            bxtr.vars[tar_name]['data_path'],
            f'{bxtr.site}_{bxtr.start_date.split("-")[0]}_{bxtr.end_date.split("-")[0]}.csv'
        )
        src_df = pd.read_csv(site_data_path)
        src_name = bxtr.vars[tar_name]['sourceVariableName']
        if src_name in src_df:
            sel_data_raw = src_df[src_name].values
            bxtr.log_var_start(tar_name)

            # convert half hourly w/m2 values to MJ/m2[hh]
            sel_data_raw = sel_data_raw * 1800 * 1e-6

            if bxtr.temporal_resolution == 'hourly':
                sel_data = sel_data_raw.reshape(-1,2).sum(1)
            elif bxtr.temporal_resolution == 'daily':
                sel_data = sel_data_raw.reshape(-1,48).sum(1)
            else:
                sys.exit(f"{bxtr.temporal_resolution} aggregation is not implemented for ONEFlux extractor in oneflux_swinpot_csvload.py")

            data = xr.Dataset({
                tar_name:
                xr.DataArray(data=sel_data,
                             dims=['time'],
                             coords={'time': date_})
            })
            data = bxtr.convert_units(data, tar_name)
            src_data.append(data)
            bxtr.log_var_end(data, tar_name, None)
        else:
            logger.warning(
                f'{bxtr.site}: target: {tar_name}, src: {src_name} was not found in the csv data frame in {site_data_path}'
            )

    src_dataset = bxtr.merge_and_format(src_data)
    return src_dataset


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